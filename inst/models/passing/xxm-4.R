# Penicillin: Random-Intercept, Cross-Classified model
# http://xxm.times.uh.edu/learn-xxm/penicillin/

libraries <- rownames(installed.packages())
if (!("lme4" %in% libraries)) stop("SKIP")
library(lme4)
library(OpenMx)

plateModel <- mxModel(
    'plateModel', type="RAM",
    latentVars = c('plate'),
    mxData(data.frame(plate=unique(Penicillin$plate)), 'raw', primaryKey='plate'),
    mxPath('plate', arrows=2, values=1))

sampleModel <- mxModel(
    'sampleModel', type="RAM",
    latentVars = c('sample'),
    mxData(data.frame(sample=unique(Penicillin$sample)), 'raw', primaryKey='sample'),
    mxPath('sample', arrows=2, values=1))

diameterModel <- mxModel(
    'diameterModel', type="RAM", plateModel, sampleModel,
    manifestVars='diameter',
    mxData(Penicillin, 'raw'),
    mxPath('one', 'diameter'),
    mxPath('diameter', arrows=2, values=1),
    mxPath('plateModel.plate', 'diameter', free=FALSE, values=1, joinKey='plate'),
    mxPath('sampleModel.sample', 'diameter', free=FALSE, values=1, joinKey='sample'))

diameterModel <- mxRun(diameterModel)

(fm03.1 <- lmer(diameter ~ 1 + (1 | plate) + (1 | sample), REML = FALSE, Penicillin))

omxCheckCloseEnough(logLik(diameterModel), logLik(fm03.1), 1e-6)
