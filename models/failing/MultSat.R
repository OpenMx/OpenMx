require(OpenMx)
require(MASS)
set.seed(200); rs=.5; xy <- mvrnorm (1000, c(0,0), matrix(c(1,rs,rs,1),2,2))
testData <- xy
summary(testData); cov(testData)
manifestVars <- c('X', 'Y')
dimnames(testData) <- list(NULL, manifestVars)

# Currently back-end is not using 'M' matrix in RAM calculation
multSatCovPatModel <- mxModel("multSatCovPat",
	manifestVars = manifestVars,
	mxPath(from="X", to="Y", arrows=2, free=T, values=.5, lbound=.01, labels="covXY"),
	mxPath(from=c("X", "Y"), arrows=2, free=T, values=1, lbound=.01, labels=c("varX","varY")),
	mxData(cov(testData), type="cov", numObs=1000, means=colMeans(testData)), type="RAM")
multSatCovPatFit <- mxRun(multSatCovPatModel)
print(multSatCovPatFit[['S']]@values); print(multSatCovPatFit@objective)

multSatRawPatModel <- mxModel("multSatRawPat",
	manifestVars = manifestVars,
	mxPath(from="X", to="Y", arrows=2, free=T, values=.5, lbound=.01, labels="covXY"),
	mxPath(from=c("X", "Y"), arrows=2, free=T, values=1, lbound=.01, labels=c("varX","varY")),
	mxData(testData, type="raw"), type="RAM")
multSatRawPatFit <- mxRun(multSatRawPatModel)
print(multSatRawPatFit[['S']]@values); print(multSatRawPatFit@objective)

#multSatCovMatModel <- mxModel("multSatCovMat",
#	mxMatrix("Symm", nrow=2, ncol=2, free=T, values=c(1,.2,.2,1), name="expCov"),
#   mxMatrix("Full", nrow=1, ncol=2, free=T, values=c(0,0), name="expMean"), 
#	mxMLObjective("expCov", "expMean"),
#	data = mxData(cov(testData), type="cov", numObs=1000, means=colMeans(testData)))
#print(multSatCovMatModel@matrices)
#multSatCovMatFit <- mxRun(multSatCovMatModel)
#multSatCovMatFit[['expCov']]@values
#print(multSatRawModel@objective)

multSatRawMatModel <- mxModel("multSatRawMat",
	mxMatrix("Symm", nrow=2, ncol=2, free=T, values=c(1,.2,.2,1), name="expCov", 
		dimnames = list(manifestVars, manifestVars)),
	mxMatrix("Full", nrow=1, ncol=2, free=T, values=c(0,0), name="expMean", 
		dimnames = list(NULL, manifestVars)),
	mxFIMLObjective("expCov", "expMean"), mxData(testData, type="raw"))
multSatRawMatFit <- mxRun(multSatRawMatModel)
print(multSatRawMatFit[['expCov']]@values); print(multSatRawMatFit[['expMean']]@values)
print(multSatRawMatFit@objective)

# Testing the last model, now with explicit parameter names on symmetric matrix
multSatRawMatModel2 <- mxModel("multSatRawMat2",
	mxMatrix("Symm", nrow=2, ncol=2, free=T, values=c(1,.2,.2,1), name="expCov", 
		labels = c('covXX', 'covXY', 'covYY'), 
		dimnames = list(manifestVars, manifestVars)),
	mxMatrix("Full", nrow=1, ncol=2, free=T, values=c(0,0), name="expMean",
		labels = c('meanX', 'meanY'),
		dimnames = list(NULL, manifestVars)),
	mxFIMLObjective("expCov", "expMean"), mxData(testData, type="raw"))
multSatRawMatFit2 <- mxRun(multSatRawMatModel2)
print(multSatRawMatFit2[['expCov']]@values); print(multSatRawMatFit2[['expMean']]@values)
print(multSatRawMatFit2@objective)

