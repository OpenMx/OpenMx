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
    # mxPath(from="one", to=manifests, arrows=1, free=T, values=colMeans(demoOneFactor)),
    # mxData(demoOneFactor, type="raw"),
    mxMatrix("Iden", nrow=nVars, name="I"),
    mxMatrix("Full", free=FALSE, values=diag(nrow=nManifest, ncol=nVars), name="Eff"),
    mxAlgebra(Eff%*%solve(I-A), name="Z"),
    mxAlgebra(Z%*%S%*%t(Z), name="C"),
    mxAlgebra(sqrt(diag2vec(C)), name="P"),
    mxCI(c("P"))
)
factorFit <- mxRun(factorModel, intervals=FALSE)
omxCheckCloseEnough(factorFit$output$fit, -3660.596, .01)

factorFitCI <- mxRun(mxModel(factorFit, mxComputeConfidenceInterval(plan=mxComputeGradientDescent())), suppressWarnings = TRUE)
factorSummCI <- summary(factorFitCI)
summary(factorFitCI)

omxCheckCloseEnough(coef(factorFit), coef(factorFitCI))
omxCheckCloseEnough(factorFit$output$fit, factorFitCI$output$fit, 0)
omxCheckCloseEnough(mxEval(Z, factorFit), mxEval(Z, factorFitCI))
omxCheckCloseEnough(mxEval(C, factorFit), mxEval(C, factorFitCI))
omxCheckCloseEnough(mxEval(P, factorFit), mxEval(P, factorFitCI))
omxCheckCloseEnough(mxEval(objective, factorFit)[1,1], mxEval(objective, factorFitCI)[1,1])

if (0) {
  options(digits=12)
  print(factorFitCI$output$fit)
  print(factorFitCI$output$computes[[2]])
}

ci <- factorFitCI$output$confidenceIntervals
#print(ci)
#cat(deparse(round(ci[,'ubound'],4)))
omxCheckCloseEnough(ci[,'estimate'], c(0.4456, 0.5401, 0.6116, 0.7302, 0.8187), .001)

lbound <- c(0.406, 0.485, 0.553, 0.6872, 0.769)
ubound <- c(0.4747, 0.5754, 0.6516, 0.778, 0.8723)

omxCheckCloseEnough(ci[,'lbound'], lbound, .03)

if (mxOption(NULL, 'Default optimizer') != "NPSOL") {
	# NPSOL needs to get slightly closer to the MLE to nail all of these
	omxCheckCloseEnough(ci[,'ubound'], ubound, .06)
}

factorParallel <- omxParallelCI(factorFit)
pci <- factorParallel$output$confidenceIntervals
omxCheckCloseEnough(pci[,'lbound'], lbound, .03)
omxCheckCloseEnough(pci[,'ubound'], ubound, .06)
