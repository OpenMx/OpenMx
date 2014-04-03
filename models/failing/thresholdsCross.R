# try to reverse the estimated thresholds by putting them very close together

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
	mxFIMLObjective("S", "M", names(oDat), "T")
	)

firstTest <- mxRun(first, unsafe=TRUE)

second <- mxModel("Sub",
	mxData(oDat, "raw"),
	mxMatrix("Stand", 2, 2, TRUE, sigma, name="S"),
	mxMatrix("Zero", 1, 2, name="M"),
	mxMatrix("Full", 2, 2, TRUE, c(0, 0.001), name="T"),
	mxFIMLObjective("S", "M", names(oDat), "T")
	)

secondTest <- mxRun(second, unsafe=TRUE)

# show that the error is due to threshold reversal
firstCross  <- any((firstTest$T$values[2,] - firstTest$T$values[1,])<0)
secondCross <- any((secondTest$T$values[2,] - secondTest$T$values[1,])<0)

# use omxCheckTrue to throw an error
omxCheckTrue(!firstCross)
omxCheckTrue(!secondCross)