library(OpenMx)

data(demoOneFactor)
factorModel <- mxModel(
    "One Factor",
    mxMatrix("Full", 5, 1, values=0.2,
	     free=TRUE, name="A"),
    mxMatrix("Symm", 1, 1, values=1,
	     free=FALSE, name="L"),
    mxMatrix("Diag", 5, 5, values=1,
	     free=TRUE, name="U"),
    mxAlgebra(A %*% L %*% t(A) + U, name="R"),
    mxExpectationNormal("R", dimnames = names(demoOneFactor)),
    mxFitFunctionML(),
    mxData(cov(demoOneFactor), type="cov", numObs=500),
    mxComputeSequence(list(
      mxComputeNumericDeriv(),
      mxComputeReportDeriv()))
    )
fitModel <- mxRun(factorModel)
fullG <- fitModel$output$gradient
fullH <- fitModel$output$hessian

limModel <- mxRun(mxModel(factorModel,
                          mxComputeSequence(list(
                            mxComputeNumericDeriv(filter=c("One Factor.A[1,1]",
                                                           "One Factor.A[2,1]")),
                            mxComputeReportDeriv()))))
omxCheckCloseEnough(limModel$output$gradient[1:2], fullG[1:2], 1e-3)
omxCheckCloseEnough(limModel$output$hessian[,1:2], fullH[,1:2], 1e-3)
