# Pastes: Three level Random-Intercepts Model
# http://xxm.times.uh.edu/learn-xxm/lme4-example-pastes/

libraries <- rownames(installed.packages())
if (!("lme4" %in% libraries)) stop("SKIP")

library(lme4)
library(OpenMx)

(fm01 <- lmer(strength ~ 1 + (1 | sample) + (1 | batch), REML = FALSE, Pastes))

batch <- mxModel(
    'batch', type="RAM",
    latentVars = c('batch'),
    mxData(data.frame(batch=unique(Pastes$batch)), 'raw', primaryKey='batch'),
    mxPath('batch', arrows=2, values=1))

sample <- mxModel(
    'sample', type="RAM", batch,
    latentVars = c('sample'),
    mxData(as.data.frame(Pastes[!duplicated(Pastes$sample), c('batch','sample')]),
           'raw', primaryKey='sample'),
    mxPath('sample', arrows=2, values=1),
    mxPath('batch.batch', 'sample', values=1, free=FALSE, joinKey='batch'))

strength <- mxModel(
    'strength', type='RAM', sample,
    manifestVars = c('strength'),
    mxData(Pastes, 'raw'),
    mxPath('one', 'strength'),
    mxPath('strength', arrows=2, values=1),
    mxPath('sample.sample', 'strength', free=FALSE, values=1, joinKey="sample"))

strength <- mxRun(strength)

omxCheckCloseEnough(strength$output$fit, 247.9945, 1e-2)

omxCheckCloseEnough(logLik(fm01), logLik(strength), 1e-6)
