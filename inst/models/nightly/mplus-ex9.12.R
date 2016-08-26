# MPlus: Two-level growth model for a continuous outcome (three-level analysis)
# https://www.statmodel.com/usersguide/chapter9.shtml

library(OpenMx)

options(width=120)
ex912 <- suppressWarnings(try(read.table("models/nightly/data/ex9.12.dat")))
if (is(ex912, "try-error")) ex912 <- read.table("data/ex9.12.dat")
colnames(ex912) <- c(paste0('y',1:4), 'x', 'w', 'clus')
ex912$clus <- as.integer(ex912$clus)

betweenModel <- mxModel(
    'betweenModel', type='RAM',
    latentVars=c(paste0('y',1:4), 'ib', 'sb', 'w'),
    mxData(ex912[!duplicated(ex912$clus),], 'raw', primaryKey='clus'),
    mxPath('ib', paste0('y',1:4), free=FALSE, values=1),
    mxPath('sb', paste0('y',1:4), free=FALSE, values=0:3),
    mxPath(c('ib','sb'), arrows=2, connect="unique.pairs", values=c(1,0,1)),
    mxPath('one', 'w', free=FALSE, labels="data.w"),
    mxPath('w', c('ib','sb')),
    mxPath('one', c('ib','sb')))

withinModel <- mxModel(
    'withinModel', type='RAM', betweenModel,
    manifestVars=paste0('y',1:4), latentVars=c('iw','sw','x'),
    mxData(ex912, 'raw'),
    mxPath('iw', paste0('y',1:4), free=FALSE, values=1),
    mxPath('sw', paste0('y',1:4), free=FALSE, values=0:3),
    mxPath(paste0('y',1:4), arrows=2, labels="yVar"),
    mxPath(c('iw','sw'), arrows=2, connect="unique.pairs", values=c(1,0,1)),
    mxPath('one', 'x', free=FALSE, labels="data.x"),
    mxPath('x', c('iw','sw')),
    mxPath(paste0('betweenModel.y', 1:4), paste0('y',1:4), free=FALSE, values=1,
	   joinKey="clus"))

withinModel <- mxRun(withinModel)

omxCheckCloseEnough(logLik(withinModel), -6531.52, 1e-2)  # matches Mplus
