require(OpenMx)

data(demoOneFactor)

manifests <- names(demoOneFactor)
latents <- c("factor")
nManifest <- length(manifests)
nVars <- nManifest + length(latents)

factorModel <- mxModel("One Factor", type="RAM",
    manifestVars = manifests,
    latentVars = latents,
    mxPath(from=latents, to=manifests, free=c(FALSE,TRUE,TRUE,TRUE,TRUE), values=1),
    mxPath(from=manifests, arrows=2, lbound=.0001),
    mxPath(from=latents, arrows=2, free=TRUE, values=1.0),
    mxData(cov(demoOneFactor), type="cov", numObs=500),
    # mxPath(from="one", to=manifests, arrows=1, free=T, values=mean(demoOneFactor)),
    # mxData(demoOneFactor, type="raw"),
    mxMatrix("Iden", nrow=nVars, name="I"),
    mxMatrix("Full", free=FALSE, values=diag(nrow=nManifest, ncol=nVars), name="Eff"),
    mxAlgebra(Eff%*%solve(I-A), name="Z"),
    mxAlgebra(Z%*%S%*%t(Z), name="C"),
    mxAlgebra(sqrt(diag2vec(C)), name="P"),
    mxCI(c("P"))
)
factorFit <- mxRun(factorModel, intervals=FALSE)
factorFitCI <- mxRun(factorFit, intervals=TRUE, suppressWarnings = TRUE)
factorSummCI <- summary(factorFitCI)
summary(factorFitCI)

omxCheckCloseEnough(factorSummCI$CI["One Factor.P[1,1]",], c(.419, .446, .474), .005)
omxCheckCloseEnough(factorSummCI$CI["One Factor.P[2,1]",], c(.508, .541, .575), .005)
omxCheckCloseEnough(factorSummCI$CI["One Factor.P[3,1]",], c(.575, .613, .651), .005)
omxCheckCloseEnough(factorSummCI$CI["One Factor.P[4,1]",], c(.687, .732, .778), .005)
omxCheckCloseEnough(factorSummCI$CI["One Factor.P[5,1]",], c(.770, .820, .872), .005)

factorParallel <- omxParallelCI(factorFit)
omxCheckCloseEnough(factorParallel@output$confidenceIntervals,
	factorFitCI@output$confidenceIntervals, 0.001)

# TODO : Compare to old Mx values.
