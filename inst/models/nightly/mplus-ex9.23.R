# MPlus: Three-level growth model with a continuous outcome and one covariate on each of the three levels
# https://www.statmodel.com/usersguide/chapter9.shtml

library(OpenMx)

#mxOption(NULL, "Number of Threads", 8L)

options(width=120)
ex923 <- suppressWarnings(try(read.table("models/nightly/data/ex9.23.dat")))
if (is(ex923, "try-error")) ex923 <- read.table("data/ex9.23.dat")
colnames(ex923) <- c(paste0('y',1:4), 'x', 'w', 'z', 'level2', 'level3')
ex923$level2 <- as.integer(ex923$level2)
ex923$level3 <- as.integer(ex923$level3)

level3Model <- mxModel(
    'level3Model', type='RAM',
    latentVars=c(paste0('y',1:4), 'ib3', 'sb3', 'z'),
    mxData(ex923[!duplicated(ex923$level3),], 'raw', primaryKey='level3'),
    mxPath('ib3', paste0('y',1:4), free=FALSE, values=1),
    mxPath('sb3', paste0('y',1:4), free=FALSE, values=0:3),
    mxPath(c('ib3','sb3'), arrows=2, connect="unique.pairs", values=c(1,0,1)),
    mxPath('one', 'z', free=FALSE, labels="data.z"),
    mxPath('z', c('ib3','sb3')),
    mxPath('one', c('ib3','sb3')))

level2Model <- mxModel(
    'level2Model', type='RAM', level3Model,
    latentVars=c(paste0('y',1:4), 'ib2', 'sb2', 'w'),
    mxData(ex923[!duplicated(ex923$level2),], 'raw', primaryKey='level2'),
    mxPath('ib2', paste0('y',1:4), free=FALSE, values=1),
    mxPath('sb2', paste0('y',1:4), free=FALSE, values=0:3),
    mxPath(c('ib2','sb2'), arrows=2, connect="unique.pairs", values=c(1,0,1)),
    mxPath('one', 'w', free=FALSE, labels="data.w"),
    mxPath('w', c('ib2','sb2')),
    mxPath(paste0('y',1:4), arrows=2),
    mxPath(paste0('level3Model.y', 1:4), paste0('y',1:4), free=FALSE, values=1,
	   joinKey="level3"))

withinModel <- mxModel(
    'withinModel', type='RAM', level2Model,
    manifestVars=paste0('y',1:4), latentVars=c('iw','sw','x'),
    mxData(ex923, 'raw'),
    mxPath('iw', paste0('y',1:4), free=FALSE, values=1),
    mxPath('sw', paste0('y',1:4), free=FALSE, values=0:3),
    mxPath(paste0('y',1:4), arrows=2),
    mxPath(c('iw','sw'), arrows=2, connect="unique.pairs", values=c(1,0,1)),
    mxPath('one', 'x', free=FALSE, labels="data.x"),
    mxPath('x', c('iw','sw')),
    mxPath(paste0('level2Model.y', 1:4), paste0('y',1:4), free=FALSE, values=1,
	   joinKey="level2"))

withinModel <- mxRun(withinModel)

omxCheckEquals(withinModel$expectation$debug$rampartUsage, c(6000, 1450))

omxCheckCloseEnough(logLik(withinModel), -56044.82, 1e-2)  # matches Mplus
