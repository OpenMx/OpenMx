library(OpenMx)

got <- suppressWarnings(try(load(file="data/heck-thomas-2015ch4.rda")))
if (is(got, "try-error")) load("models/passing/data/heck-thomas-2015ch4.rda")

orgModel <- mxModel(
  "org", type="RAM",
  mxData(level2, 'raw', primaryKey = "orgid"),
  manifestVars = c('zproduct'),
  latentVars = c('Benefits','Conditions', 'org1','org2'),
  mxPath(c('zproduct','Benefits','Conditions'), arrows=2, values=1),
  mxPath('Benefits','Conditions', arrows=2),
  mxPath('one', 'org1', labels='data.org1', free=FALSE),
  mxPath('one', 'org2', labels='data.org2', free=FALSE),
  mxPath('org1', c('Benefits', 'Conditions')),
  mxPath('org2', 'zproduct', labels="a"),
  mxPath('org2', c('Benefits', 'Conditions')),
  mxPath('one', 'zproduct'),
  mxPath('zproduct', c('Benefits','Conditions'), labels=paste0('b',1:2)),
  mxMatrix(nrow=2, ncol=1, free=TRUE, labels=paste0('b',1:2), name='B'),
  mxAlgebra(B * a, name="indirect")
)

empModel <- mxModel(
  "emp", type="RAM", orgModel,
  mxData(level1, 'raw'),
  manifestVars = c('benefit','cond'),
  latentVars = c('female','white'),
  mxPath('one','female', labels='data.female', free=FALSE),
  mxPath('one','white', labels='data.white', free=FALSE),
  mxPath(c('female','white'), c('benefit','cond'), connect = "all.bivariate"),
  mxPath(c('benefit','cond'), arrows=2, connect = "unique.pairs", values=c(1,0,1)),
  mxPath('org.Benefits', 'benefit', values=1, free=FALSE, joinKey = "orgid"),
  mxPath('org.Conditions', 'cond', values=1, free=FALSE, joinKey = "orgid"),
  mxPath('one', c('benefit','cond')),
  mxCI("org.indirect")
)

empModel <- omxSetParameters(empModel, names(coef(empModel)),
                             values=c(0.06, 0.06, -0.01, -0.01, 2.16, 1.48, 2.09, 4.82,
                                      5.1, 0.33, 0.31, 0.29, 0.31, -0.45, 0.05, 0.04,
                                      1.07, 0.13, 0.09,  0.09, 0.09))

empModel <- mxRun(mxModel(empModel,
                          mxComputeOnce('fitfunction', 'fit')))

empModel$expectation$.useSufficientSets <- FALSE
empModel2 <- mxRun(mxModel(empModel,
                           mxComputeOnce('fitfunction', 'fit')))

omxCheckCloseEnough(empModel$output$fit, empModel2$output$fit, .01)
