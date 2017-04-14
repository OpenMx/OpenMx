# ===========
# = history =
# ===========
# 2017-04-14 05:56PM TBATES
# update with 	mxExpectationNormal, mxFitFunctionML()

# Error is that OpenMx allows estimated thresholds to get out of order.
# This script trigger this by starting them very close together

library(OpenMx)
library(mvtnorm)

set.seed(2)

n     <- 10000
nvar  <- 2
k     <- .5
sigma <- matrix(k, nvar, nvar)
diag(sigma) <- 1
tloc  <- round(n/2)


# generate data
cDat  <- rmvnorm(n, sigma=sigma)
oDat  <- (cDat>0) + (cDat>10/n)

# make sure that no category is empty
any(oDat[,1]==1)&any(oDat[,2]==1)

# mxFactor the data
oDat  <- mxFactor(data.frame(oDat), 0:2)

# take a look
summary(oDat)

# cause the error!
first <- mxModel("Sub",
	mxData(oDat, "raw"),
	mxMatrix("Stand", 2, 2, TRUE, sigma, name="S"),
	mxMatrix("Zero", 1, 2, name="M"),
	mxMatrix("Full", 2, 2, TRUE, c(-1, 1), name="T"),
	mxExpectationNormal("S", "M", names(oDat), "T"),
	mxFitFunctionML()
)

firstTest <- mxRun(first, unsafe=TRUE)
# as of 2.7.9 gives code red
# In model 'Sub' Optimizer returned a non-zero status code 6. The model does not satisfy the first-order optimality conditions to the required accuracy, and no improved point for the merit function could be found during the final linesearch (Mx status RED)

second <- mxModel("Sub",
	mxData(oDat, "raw"),
	mxMatrix("Stand", 2, 2, TRUE, sigma, name="S"),
	mxMatrix("Zero", 1, 2, name="M"),
	mxMatrix("Full", 2, 2, TRUE, c(0, 0.001), name="T"),
	mxExpectationNormal("S", "M", names(oDat), "T"),
	mxFitFunctionML()
)

secondTest <- mxRun(second, unsafe=TRUE)
# Warning messages:
# 1: In model 'Sub' Optimizer returned a non-zero status code 10. Starting values are not feasible. Consider mxTryHard()
# 2: The job for model 'Sub' exited abnormally with the error message: fit is not finite (0: Found 1 thresholds too close together in column 1.
# 1: Found 1 thresholds too close together in column 1.
# 2: Found 1 thresholds too close together in column 1.
# 3: Found 1 thresholds too close together in column 1.
# )

# show that the error is due to threshold reversal
firstCross  <- any((firstTest$T$values[2,] - firstTest$T$values[1,])<0)
secondCross <- any((secondTest$T$values[2,] - secondTest$T$values[1,])<0)

# use omxCheckTrue to throw an error
omxCheckTrue(!firstCross)
omxCheckTrue(!secondCross)