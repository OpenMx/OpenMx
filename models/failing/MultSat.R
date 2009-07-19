require(OpenMx)
require(MASS)
set.seed(200); rs=.5; xy <- mvrnorm (1000, c(0,0), matrix(c(1,rs,rs,1),2,2))
testData <- xy
summary(testData); cov(testData)
dimnames(testData) <- list(NULL, c('X','Y'))

# Currently back-end is not using 'M' matrix in RAM calculation
multSatCovPatModel <- mxModel("multSatCovPat",
	manifestVars = c("X", "Y"),
	mxPath(from="X", to="Y", arrows=2, free=T, values=.5, lbound=.01, labels="covXY"),
	mxPath(from=c("X", "Y"), arrows=2, free=T, values=1, lbound=.01, labels=c("varX","varY")),
	mxData(cov(testData), type="cov", numObs=1000, means=colMeans(testData)), type="RAM")
multSatCovPatFit <- mxRun(multSatCovPatModel)
print(multSatCovPatFit[['S']]@values); print(multSatCovPatFit@objective)

multSatRawPatModel <- mxModel("multSatRawPat",
	manifestVars = c("X", "Y"),
	mxPath(from="X", to="Y", arrows=2, free=T, values=.5, lbound=.01, labels="covXY"),
	mxPath(from=c("X", "Y"), arrows=2, free=T, values=1, lbound=.01, labels=c("varX","varY")),
	mxData(testData, type="raw"), type="RAM")
multSatRawPatFit <- mxRun(multSatRawPatModel)
print(multSatRawPatFit[['S']]@values); print(multSatRawPatFit@objective)

# Crashes omxRAMObjective.c
#multSatCovMatModel <- mxModel("multSatCovMat",
#	mxMatrix("Symm", nrow=2, ncol=2, free=T, values=c(1,.2,.2,1), name="expCov"),
#	mxRAMObjective(),
#	data = mxData(cov(testData), type="cov", numObs=1000, means=colMeans(testData)), type="RAM")
#print(multSatCovMatModel@matrices)
#multSatCovMatFit <- mxRun(multSatCovMatModel)
#multSatCovMatFit[['expCov']]@values
#print(multSatRawModel@objective)

multSatRawMatModel <- mxModel("multSatRawMat",
	mxMatrix("Symm", nrow=2, ncol=2, free=T, values=c(1,.2,.2,1), name="expCov", 
		dimnames = list(c('X','Y'), c('X','Y'))),
	mxMatrix("Full", nrow=1, ncol=2, free=T, values=c(0,0), name="expMean", 
		dimnames = list(NULL, c('X', 'Y'))),
	mxFIMLObjective("expCov", "expMean"), mxData(testData, type="raw"))
multSatRawMatFit <- mxRun(multSatRawMatModel)
print(multSatRawMatFit[['expCov']]@values); print(multSatRawMatFit[['expMean']]@values)
print(multSatRawMatFit@objective)
