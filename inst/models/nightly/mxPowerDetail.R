library(OpenMx)

data(demoOneFactor)
manifests <- names(demoOneFactor)
latents <- c("G")
factorModel <- mxModel("One Factor",
                       type="RAM",
                       manifestVars = manifests,
                       latentVars = latents,
                       mxPath(from=latents, to=manifests, values=0.8),
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

got <- mxPower(factorModelFit, indModel, power=NA, n=1/5, method='ncp')
detail <- attr(got, 'detail')
omxCheckEquals(sort(names(detail)),
               c("method", "n", "power", "sig.level", "statistic"))
omxCheckEquals(detail$method, 'ncp')
omxCheckEquals(detail$n, 1/5)
omxCheckCloseEnough(detail$power, .9, .01)
omxCheckEquals(detail$sig.level, .05)
omxCheckEquals(detail$statistic, 'LRT')

got <- mxPower(factorModelFit, indModel, power=.8, method='ncp')
detail <- attr(got, 'detail')
omxCheckEquals(sort(names(detail)),
               c("method", "n", "power", "sig.level", "statistic"))
omxCheckEquals(detail$method, 'ncp')
omxCheckEquals(detail$n, .159, 1e-3)
omxCheckEquals(detail$power, .8)
omxCheckEquals(detail$sig.level, .05)
omxCheckEquals(detail$statistic, 'LRT')

got <- mxPower(factorModelFit, indModel, probes=30, power=.8, n=.08)
detail <- attr(got, 'detail')
omxCheckEquals(sort(names(detail)),
               c("method", "n", "parameter", "parameterDiff", "power", "probes",  "sig.level", "statistic"))
omxCheckEquals(detail$method, 'empirical')
omxCheckEquals(detail$n, .08)
omxCheckEquals(detail$parameter, 'One Factor.A[1,6]')
omxCheckCloseEnough(detail$parameterDiff, 0.154, .03)
omxCheckEquals(detail$power, .8)
omxCheckEquals(detail$probes, 30)
omxCheckEquals(detail$sig.level, .05)
omxCheckEquals(detail$statistic, 'LRT')
