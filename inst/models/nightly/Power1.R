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
indModel <- mxRun(indModel)

set.seed(1)

got <- mxPower(factorModelFit, indModel)

got <- mxPower(factorModelFit, indModel, prevOutput = got,
               grid=seq(15,160,length.out = 20))

got2 <- mxPower(factorModelFit, indModel, method='ncp',
                grid=seq(15,160,length.out = 20))

omxCheckCloseEnough(c(pmin(got2[,'p'] - got[,'pmin'], 0),
                      pmin(got[,'pmax'] - got2[,'p'], 0)),
                    rep(0,40), .05)

# --------------------

indModel <- factorModelFit
indModel$A$values['x1','G'] <- .1
indModel$A$free['x1','G'] <- FALSE
indModel <- mxRun(indModel)

got <- mxPower(factorModelFit, indModel, probes = 50, n=100)
got <- mxPower(factorModelFit, indModel, n=100, prevOutput = got)
omxCheckCloseEnough(got[findInterval(.8, got$p), 'loading1'],
                    .1285, .005)
