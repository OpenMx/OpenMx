# stuff to dye for
# http://xxm.times.uh.edu/learn-xxm/lme4-example-dyestuff/

libraries <- rownames(installed.packages())
if (!("lme4" %in% libraries)) stop("SKIP")

library(lme4)
library(OpenMx)

batch <- mxModel(
    'batch', type="RAM",
    latentVars = c('batch'),
    mxData(data.frame(batch=unique(Dyestuff$Batch)), 'raw', primaryKey='batch'),
    mxPath('batch', arrows=2, values=100, ubound=10000))

yield <- mxModel(
    'yield', type='RAM', batch,
    manifestVars = c('Yield'),
    mxData(Dyestuff, 'raw'),
    mxPath('one', 'Yield'),
    mxPath('Yield', arrows=2, values=100, ubound=10000),
    mxPath('batch.batch', 'Yield', free=FALSE, values=1, joinKey="Batch"))

yield <- mxRun(yield)

omxCheckCloseEnough(yield$output$fit, 327.3271, 1e-2)

(fm01 <- lmer(Yield ~ 1 + (1 | Batch), REML = FALSE, Dyestuff))
omxCheckCloseEnough(logLik(fm01), logLik(yield), 5e-3)

# Check free parameter estimates
lVars <- as.data.frame(VarCorr(fm01))[, 'vcov']
lEst <- c(lVars[2], fixef(fm01), lVars[1])
oEst <- coef(yield)
rms <- function(x, y){sqrt(mean((x-y)^2))}
omxCheckCloseEnough(rms(oEst, lEst), 0, 3)
