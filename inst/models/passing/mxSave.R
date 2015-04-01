library(OpenMx)

# Simulate some data

sampleSize <- 250
x=rnorm(sampleSize, mean=0, sd=1)
y= 0.5*x + rnorm(sampleSize, mean=0, sd=1)
tmpFrame <- data.frame(x, y)
tmpNames <- names(tmpFrame)

# Create a model that includes an expected covariance matrix,
# an expectation function, a fit function, and an observed covariance matrix

data <- mxData(cov(tmpFrame), type="cov", numObs = sampleSize)
expCov <- mxMatrix(type="Symm", nrow=2, ncol=2, values=c(.2,.1,.2), free=TRUE, name="expCov")
expFunction <- mxExpectationNormal(covariance="expCov", dimnames=tmpNames)
fitFunction <- mxFitFunctionML()
testModel <- mxModel(model="testModel", expCov, data, expFunction, fitFunction)

#Use mxRun to optimize the free parameters in the expected covariance matrix
modelOut <- mxRun(testModel, checkpoint = TRUE)
modelOut$expCov

# Save the ending state of modelOut in a checkpoint file
mxSave(modelOut)

modelRestored <- mxRestore(testModel)
omxCheckCloseEnough(modelRestored$expCov$values, modelOut$expCov$values, 1e-5)
