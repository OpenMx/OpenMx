library(OpenMx)

data(demoOneFactor)
manifests <- names(demoOneFactor)
latents <- c("G")
factorModel <- mxModel("One Factor",
                       type="RAM",
                       manifestVars = manifests,
                       latentVars = latents,
                       mxPath(from=latents, to=manifests, values=0.8,
                              labels=paste0("loading", 1:length(manifests))),
                       mxPath(from=manifests, arrows=2,values=1),
                       mxPath(from=latents, arrows=2,
                              free=FALSE, values=1.0),
                       mxPath(from="one", to=manifests),
                       mxData(demoOneFactor, type="raw"))
factorModelFit <- mxRun(factorModel)

indModel <- factorModelFit
indModel$A$values['x1','G'] <- 0.3
indModel$A$free['x1','G'] <- FALSE

got4 <- mxPowerSearch(factorModelFit, indModel, method = 'ncp')
omxCheckCloseEnough(got4[findInterval(.8, got4$power), 'N'],
                    34.67, 1)

indModel <- mxRun(indModel)

set.seed(1)

got <- mxPowerSearch(factorModelFit, indModel)

got <- mxPowerSearch(factorModelFit, indModel, previousRun = got,
               grid=seq(15,160,length.out = 20))

got2 <- mxPowerSearch(factorModelFit, indModel, method='ncp',
                grid=seq(15,160,length.out = 20))

omxCheckCloseEnough(c(pmin(got2[,'power'] - got[,'lower'], 0),
                      pmin(got[,'upper'] - got2[,'power'], 0)),
                    rep(0,40), .05)

# --------------------

refs <- mxRefModels(factorModelFit, run = TRUE)
mxCompare(refs[['Saturated']], factorModelFit)
got3 <- mxPowerSearch(factorModelFit, refs[['Saturated']],
                statistic = 'AIC', probes = 300)
omxCheckCloseEnough(got3[findInterval(.8, got3$power), 'N'],
                    16.45, 3)

# --------------------

indModel <- factorModelFit
indModel$A$values['x1','G'] <- .1
indModel$A$free['x1','G'] <- FALSE
indModel <- mxRun(indModel)

got <- mxPowerSearch(factorModelFit, indModel, probes = 50, n=100)
got <- mxPowerSearch(factorModelFit, indModel, n=100, previousRun = got)
omxCheckCloseEnough(got[findInterval(.8, got$power), 'loading1'],
                    .1285, .005)
