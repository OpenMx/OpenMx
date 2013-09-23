
require(OpenMx)
data(demoOneFactor)
demoOneFactorMiss <- as.matrix(demoOneFactor)[1:10,]
nvar <- ncol(demoOneFactor)
varnames <- colnames(demoOneFactor)



missmat <- matrix(c(2, 4, 5, 1, 2, 3), nrow=3, ncol=2)
demoOneFactorMiss[missmat] <- NA

ssModel <- mxModel(name="State Space Missing Debug",
	mxMatrix("Full", 1, 1, FALSE, .3, name="A"),
	mxMatrix("Zero", 1, 1, name="B"),
	mxMatrix("Full", nvar, 1, FALSE, c(.4, .5, .6, .7, .8), name="C", dimnames=list(varnames, "F1")),
	mxMatrix("Zero", nvar, 1, name="D"),
	mxMatrix("Diag", 1, 1, FALSE, 1, name="Q"),
	mxMatrix("Diag", nvar, nvar, FALSE, .2, name="R"),
	mxMatrix("Zero", 1, 1, name="x0"),
	mxMatrix("Diag", 1, 1, FALSE, 1, name="P0"),
	mxData(observed=demoOneFactorMiss, type="raw"),
	imxExpectationStateSpace("A", "B", "C", "D", "Q", "R", "x0", "P0"),
	mxFitFunctionML()
)

ssRun <- mxRun(ssModel)

