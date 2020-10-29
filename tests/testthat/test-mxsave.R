#
#   Copyright 2007-2019 by the individuals mentioned in the source code history
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
library(testthat)
context("mxSave")

dir <- tempdir()  # safe place to create files
mxOption(key="Checkpoint Directory", value=dir)

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
plan <- omxDefaultComputePlan()
plan$steps <- c(plan$steps,
                CK=mxComputeCheckpoint(path=file.path(dir,"backendChkpt.omx"),
                                    standardErrors = TRUE, gradient = TRUE, vcov=TRUE,
                                    toReturn=TRUE))
testModel <- mxModel(model="testFish", expCov, data,
                     mxExpectationNormal(covariance="expCov", dimnames=tmpNames),
                     mxFitFunctionML(),
                     plan)

#Use mxRun to optimize the free parameters in the expected covariance matrix
modelOut <- mxRun(testModel, checkpoint = TRUE)

modelRestored <- mxRestore(testModel, strict = TRUE)  # old checkpoint output
omxCheckCloseEnough(coef(modelRestored), coef(modelOut), 1e-5)

bckpt <- read.table(file.path(dir,"backendChkpt.omx"),
           header=TRUE, sep="\t", stringsAsFactors=FALSE, check.names=FALSE)
bckpt2 <- modelOut$compute$steps$CK$log
omxCheckEquals(colnames(bckpt), colnames(bckpt2))
for (cx in colnames(bckpt)) {
  if (cx == 'timestamp') next
  omxCheckEquals(substr(as.character(bckpt[1,cx]), 1,10),
                 substr(as.character(bckpt2[1,cx]), 1,10))
}

modelRestored <- mxRestoreFromDataFrame(testModel, bckpt)
omxCheckCloseEnough(coef(modelRestored), coef(modelOut), 1e-5)
omxCheckCloseEnough(modelRestored$output$standardErrors, modelOut$output$standardErrors, 1e-5)
omxCheckCloseEnough(modelRestored$output$gradient, modelOut$output$gradient, 1e-5)
omxCheckCloseEnough(modelRestored$output$vcov, modelOut$output$vcov, 1e-5)
summary(modelRestored)

# Save the ending state of modelOut in a checkpoint file
mxSave(modelOut, chkpt.prefix="z")

modelRestored <- mxRestore(testModel, chkpt.prefix="z", strict=TRUE)
omxCheckCloseEnough(coef(modelRestored), coef(modelOut), 1e-5)

