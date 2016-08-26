# http://xxm.times.uh.edu/learn-xxm/two-level-confirmatory-factor-analysis-with-a-random-slope/

library(OpenMx)

options(width=120)
got <- suppressWarnings(try(load("models/nightly/data/lranslp.xxm.RData")))
if (is(got, "try-error")) load("data/lranslp.xxm.RData")

for (col in paste0('l', 1:2)) {
	l1[[col]] <- as.integer(l1[[col]])
}

l2Model <- mxModel(
    "l2Model", type="RAM",
    latentVars=paste0('eta', 1:2),
    mxData(l1[!duplicated(l1$l2),], 'raw', primaryKey='l2'),
    mxPath('one', 'eta2'),
    mxPath(paste0('eta',1:2), arrows=2, connect='unique.pairs', values=c(1,0,1)))
    
l1Model <- mxModel(
    "l1Model", type="RAM", l2Model,
    manifestVars=paste0('y',1:4), latentVars='eta',
    mxData(l1, 'raw'),
    mxPath('one', paste0('y',1:4)),
    mxPath(c('eta', paste0('y',1:4)), arrows=2, values=1),
    mxPath('eta', paste0('y',1:4), free=c(FALSE, rep(TRUE,3)), values=1),
    mxPath(paste0('l2Model.eta',1:2), 'eta', free=FALSE,
	   values=c(1,0), labels=c(NA,"data.x"), joinKey='l2'))

l1Model <- mxRun(l1Model)

omxCheckCloseEnough(l1Model$output$fit, 2473.621, 1e-2)
