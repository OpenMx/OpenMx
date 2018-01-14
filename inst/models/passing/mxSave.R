#
#   Copyright 2007-2018 The OpenMx Project
#
#   Licensed under the Apache License, Version 2.0 (the "License");
#   you may not use this file except in compliance with the License.
#   You may obtain a copy of the License at
# 
#        http://www.apache.org/licenses/LICENSE-2.0
# 
#   Unless required by applicable law or agreed to in writing, software
#   distributed under the License is distributed on an "AS IS" BASIS,
#   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#   See the License for the specific language governing permissions and
#   limitations under the License.


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
