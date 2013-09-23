
require(OpenMx)
data(demoOneFactor)
demoOneFactorMiss <- as.matrix(demoOneFactor)
nvar <- ncol(demoOneFactor)
varnames <- colnames(demoOneFactor)
nmiss <- 20

set.seed(20)

missmat <- matrix(c(
	sample(1:nrow(demoOneFactor), size=nmiss),
	sample(1:nvar, size=nmiss, replace=TRUE)), ncol=2, nrow=nmiss)
demoOneFactorMiss[missmat] <- NA

ssModel <- mxModel(name="State Space Manual Example",
	mxMatrix("Full", 1, 1, TRUE, .3, name="A"),
	mxMatrix("Zero", 1, 1, name="B"),
	mxMatrix("Full", nvar, 1, TRUE, .6, name="C", dimnames=list(varnames, "F1")),
	mxMatrix("Zero", nvar, 1, name="D"),
	mxMatrix("Diag", 1, 1, FALSE, 1, name="Q"),
	mxMatrix("Diag", nvar, nvar, TRUE, .2, name="R"),
	mxMatrix("Zero", 1, 1, name="x0"),
	mxMatrix("Diag", 1, 1, FALSE, 1, name="P0"),
	mxData(observed=demoOneFactor, type="raw"),
	imxExpectationStateSpace("A", "B", "C", "D", "Q", "R", "x0", "P0"),
	mxFitFunctionML()
)
ssRun <- mxRun(ssModel)
summary(ssRun)

ssMiss <- mxModel(ssModel, name="With Missing",
	mxData(observed=demoOneFactorMiss, type="raw")
	)

ssMissRun <- mxRun(ssMiss)

ssFactor <- mxModel(ssModel, name="As Factor",
	mxMatrix("Full", 1, 1, FALSE, 0, name="A")
	)

ssFactorRun <- mxRun(ssFactor)

liFactor <- mxModel(ssFactor, name="LISREL Factor",
	mxMatrix("Full", 1, 1, FALSE, 0, name="KA", dimnames=list("F1", "F1")),
	mxMatrix("Full", nvar, 1, TRUE, 0, name="TX", dimnames=list(varnames, NA)),
	mxExpectationLISREL(LX="C", PH="Q", TD="R", TX="TX", KA="KA")
	)
liRun <- mxRun(liFactor)
summary(liRun)


summary(ssFactorRun)$parameters[, c(5, 6)]
summary(liRun)$parameters[, c(5, 6)]


ssFactorMiss <- mxModel(ssFactor, name="As Factor with Missing",
	mxData(observed=demoOneFactorMiss, type="raw")
	)

liFactorMiss <- mxModel(liFactor, name="LISREL Factor with Missing",
	mxData(observed=demoOneFactorMiss, type="raw")
	)

ssFactorMissRun <- mxRun(ssFactorMiss)
liMissRun <- mxRun(liFactorMiss)


summary(ssFactorMissRun)$parameters[, c(5, 6)]
summary(liMissRun)$parameters[, c(5, 6)]



