library(OpenMx)

# sample size and threshold
n <- 500
thresh <- 0

# make some data
x <- rnorm(n)
z <- .8 * x + rnorm(n, 0, .6)
y <- rep(0, n)
y[z>thresh] <- 1
y <- as.ordered(y)
testData <- data.frame(x, y)

# make solve(I-A) = I + A matrices
iaMat1 <- mxMatrix("Full", 2, 2,
		free=c(FALSE, FALSE, TRUE, FALSE),
		values=c(1, 0, 0, 1),
		labels=c(NA, NA, "simReg", NA),
		name="IA")
iaMat2 <- mxMatrix("Full", 2, 2,
		free=c(FALSE, TRUE, FALSE, FALSE),
		values=c(1, 0, 0, 1),
		labels=c(NA, "logReg", NA, NA),
		name="IA")

# make a model with just a regression of continuous on binary
testModel <- mxModel("JointRegressionTest",
	mxData(testData, "raw"),
	iaMat1,
	mxMatrix("Diag", 2, 2,
		free=c(TRUE, FALSE),
		values=1,
		labels=c("varx", "residy"),
		name="S"),
	mxAlgebra(IA %*% S %*% t(IA), name="cov"),
	mxMatrix("Full", 1, 2, free=c(TRUE, FALSE), 
		labels=c("meanX", NA), name="mean"),
	mxMatrix("Full", 1, 1, free=TRUE, 
		labels="threshY", dimnames=list(NA, "y"), 
		name="thresh"),
	mxFIMLObjective("cov", "mean", names(testData), thresholds="thresh", threshnames="y")
	)
	
# run it
testResults <- mxRun(testModel)

# uh-oh
summary(testResults)

# check out the Hessian
testResults@output

# here's the simple regression, which may have a scaling difference
summary(lm(x~y))