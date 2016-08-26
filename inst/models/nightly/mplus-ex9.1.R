# MPLUS: TWO-LEVEL REGRESSION ANALYSIS FOR A CONTINUOUS DEPENDENT VARIABLE WITH A RANDOM INTERCEPT
# https://www.statmodel.com/usersguide/chapter9.shtml

library(OpenMx)

options(width=120)
ex91 <- suppressWarnings(try(read.table("models/nightly/data/ex9.1a.dat")))
if (is(ex91, "try-error")) ex91 <- read.table("data/ex9.1a.dat")
colnames(ex91) <- c('y', 'x', 'w', 'xm', 'clus')
ex91$clus <- as.integer(ex91$clus)

betweenModel <- mxModel(
    'betweenModel', type='RAM',
    latentVars=c('w', 'xm', 'y'),
    mxData(ex91[!duplicated(ex91$clus),], 'raw', primaryKey='clus'),
    mxPath('y', arrows=2),
    mxPath('one', 'w', free=FALSE, labels='data.w'),
    mxPath('one', 'xm', free=FALSE, labels='data.xm'),
    mxPath(c('w', 'xm'), 'y'))
	
withinModel <- mxModel(
    'withinModel', type='RAM', betweenModel,
    mxData(ex91, 'raw'),
    manifestVars='y', latentVars=c('x'),
    mxPath('one', 'x', free=FALSE, labels='data.x'),
    mxPath('x', 'y'),
    mxPath('one', 'y'),
    mxPath('y', arrows=2, values=1),
    mxPath('betweenModel.y', 'y', free=FALSE, values=1, joinKey='clus'))

withinModel <- mxRun(withinModel)

omxCheckEquals(withinModel$expectation$debug$rampartUsage, 890)

omxCheckCloseEnough(summary(withinModel)$informationCriteria['AIC:','par'], 3063.876, 1e-3)
