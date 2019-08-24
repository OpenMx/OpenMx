library(OpenMx)
library(testthat)

suppressWarnings(RNGversion("3.5"))
set.seed(1)

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

p1 <- mxPower(factorModelFit, indModel, method = 'ncp', sig.level = .01)
expect_equivalent(c(p1), 118)
p2 <- mxPower(factorModelFit, indModel, method = 'ncp', sig.level = .005)
expect_equivalent(c(p2), 132)
expect_error(mxPower(factorModelFit, indModel, method = 'ncp',
                     sig.level = .5),
             "Try method='empirical'")

got4 <- mxPowerSearch(factorModelFit, indModel, method = 'ncp')
omxCheckCloseEnough(got4[findInterval(.8, got4$power), 'N'], 73)

# only 1 p-value handled at a time
expect_error(mxPower(factorModelFit, indModel, method = 'ncp', sig.level=c(.05,.01)),
	     "one sig.level at a time")

# more than 1 power request is fine
p3 <- mxPower(factorModelFit, indModel, method = 'ncp', sig.level = .05,power=c(.5,.8))
expect_equivalent(p3, c(44,80))

set.seed(1)

got <- mxPowerSearch(factorModelFit, indModel)
m1 <- attr(got,'model')

got <- mxPowerSearch(factorModelFit, indModel, previousRun = got,
               grid=seq(15,150,length.out = 20))
m2 <- attr(got,'model')
expect_equal(coef(m1), coef(m2), scale=1, tolerance=1e-14)

got2 <- mxPowerSearch(factorModelFit, indModel, method='ncp',
                grid=seq(15,150,length.out = 20))

omxCheckCloseEnough(c(pmin(got2[,'power'] - got[,'lower'], 0),
                      pmin(got[,'upper'] - got2[,'power'], 0)),
                    rep(0,40), .05)

# --------------------

indModel <- factorModelFit
indModel$A$values['x1','G'] <- .1
indModel$A$free['x1','G'] <- FALSE
indModel <- mxRun(indModel)

got <- mxPowerSearch(factorModelFit, indModel, probes = 50, n=100)
got <- mxPowerSearch(factorModelFit, indModel, n=100, previousRun = got)
omxCheckCloseEnough(got[findInterval(.8, got$power), 'loading1'],
                    .16, .01)

