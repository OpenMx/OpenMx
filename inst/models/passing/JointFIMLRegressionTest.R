# test of the joint ordinal/continuous fiml optimizer
# create a single continuous and single ordinal variable
# test both logistic regression and simple regression

library(OpenMx)
set.seed(1110221056)

# sample size and threshold
n <- 500
thresh <- 0

# make some data
x <- rnorm(n)
z <- .8 * x + rnorm(n, 0, .6)
y <- rep(0, n)
y[z>thresh] <- 1
#y <- as.ordered(y)
testData <- data.frame(x, as.ordered(y))
names(testData) <- c("x", "y")

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
regModel1 <- mxModel("JointSimpleRegressionTest",
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
	mxFitFunctionML(),mxExpectationNormal("cov", "mean", dimnames = names(testData), thresholds="thresh", threshnames="y")
	)
# now a logistic regression
regModel2 <- mxModel(regModel1, iaMat2, name="JointLogisticRegressionTest")
	
# run them
regResults1 <- mxRun(regModel1)
regResults2 <- mxRun(regModel2)

# the summaries
summary(regResults1)
summary(regResults2)

# check the likelihoods
omxCheckCloseEnough(regResults1$output$Minus2LogLikelihood, 1789.092, 0.01)
omxCheckCloseEnough(regResults2$output$Minus2LogLikelihood, 1789.092, 0.01)

# check the parameters
omxCheckCloseEnough(regResults1$output$estimate, c(0.847, 0.282, -0.036, 0.062), 0.001)
omxCheckCloseEnough(regResults2$output$estimate, c(1.592, 1.001, -0.036, 0.117), 0.01)
