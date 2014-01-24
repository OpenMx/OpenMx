require(OpenMx)

data(demoOneFactor)

manifests <- names(demoOneFactor)
latents <- c("factor")
nVar <- length(manifests)
nFac <- length(latents)

factorModel <- mxModel("One Factor ML",
    mxData(cov(demoOneFactor), type="cov", numObs=500),
    # mxData(demoOneFactor, type="raw"),
    mxMatrix("Full", 1, nVar, free=T, values=colMeans(demoOneFactor), name="M"),
    mxMatrix("Full", nVar, nFac, free=T, values=.2, name="A"),
    mxMatrix("Diag", nVar, nVar, free=T, values=1, lbound=.0001, name="D"),
    mxAlgebra(A%*%t(A) + D, name="C"),
    mxAlgebra(sqrt(diag2vec(C)),name="P"),
    mxMLObjective("C", dimnames=manifests),
    mxCI(c("P"))
)
factorFitCI <- mxRun(factorModel, intervals=TRUE, suppressWarnings = TRUE)
factorSummCI <- summary(factorFitCI)

factorModelRaw <- mxModel(factorFitCI, mxData(demoOneFactor, type="raw"), mxFIMLObjective("C", "M", dimnames=manifests), name = "One Factor FIML")
factorFitRawCI <- mxRun(factorModelRaw, intervals=TRUE, suppressWarnings = TRUE)
factorSummRawCI <- summary(factorFitRawCI)

omxCheckCloseEnough(factorSummCI$CI["One Factor ML.P[1,1]",], c(.419, .446, .474), .005)
omxCheckCloseEnough(factorSummCI$CI["One Factor ML.P[2,1]",], c(.508, .541, .575), .005)
omxCheckCloseEnough(factorSummCI$CI["One Factor ML.P[3,1]",], c(.575, .613, .651), .005)
omxCheckCloseEnough(factorSummCI$CI["One Factor ML.P[4,1]",], c(.687, .732, .778), .005)
omxCheckCloseEnough(factorSummCI$CI["One Factor ML.P[5,1]",], c(.770, .820, .872), .005)


omxCheckCloseEnough(factorSummRawCI$CI["One Factor FIML.P[1,1]",], c(.419, .446, .474), .005)
omxCheckCloseEnough(factorSummRawCI$CI["One Factor FIML.P[2,1]",], c(.508, .541, .575), .005)
omxCheckCloseEnough(factorSummRawCI$CI["One Factor FIML.P[3,1]",], c(.575, .613, .651), .005)
omxCheckCloseEnough(factorSummRawCI$CI["One Factor FIML.P[4,1]",], c(.687, .732, .778), .005)
omxCheckCloseEnough(factorSummRawCI$CI["One Factor FIML.P[5,1]",], c(.770, .820, .872), .005)

# Compare to original MX Estimates
#          5  Confidence intervals requested in group            1
# Matrix Element Int.      Estimate         Lower         Upper  Lfail Ufail
# P   1   1   1  95.0         0.4465        0.4192        0.4749 0 1   0 1
# P   1   1   2  95.0         0.5412        0.5081        0.5756 0 1   6 1
# P   1   1   3  95.0         0.6131        0.5753        0.6518 0 1   4 1
# P   1   1   4  95.0         0.7322        0.6870        0.7783 6 1   0 0
# P   1   1   5  95.0         0.8209        0.7702        0.8726 0 1   0 0
