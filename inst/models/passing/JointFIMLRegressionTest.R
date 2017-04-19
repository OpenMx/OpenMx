#
#   Copyright 2007-2017 The OpenMx Project
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
	mxFitFunctionML(jointConditionOn='continuous'),
	mxExpectationNormal("cov", "mean", dimnames = names(testData), thresholds="thresh", threshnames="y")
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
condOnContinuous <- 1789.092
condOnOrdinal <- 1791.779
omxCheckCloseEnough(regResults1$output$Minus2LogLikelihood, condOnContinuous, 0.01)
omxCheckCloseEnough(regResults2$output$Minus2LogLikelihood, condOnContinuous, 0.01)

# check the parameters
omxCheckCloseEnough(regResults1$output$estimate, c(0.847, 0.282, -0.036, 0.062), 0.001)
omxCheckCloseEnough(regResults2$output$estimate, c(1.592, 1.001, -0.036, 0.117), 0.01)

# Simulate similar data
test2Data <- mxGenerateData(testData, 500)
reg2Results1 <- mxRun(mxModel(regModel1, mxData(test2Data, "raw")))
sum2 <- fivenum(abs(reg2Results1$output$estimate - regResults1$output$estimate))
omxCheckTrue(sum2[1] > 0)
omxCheckTrue(sum2[5] < .1)

reg2Results1 <- mxBootstrap(reg2Results1, 10,
                            OK=c("OK", "OK/green", "nonzero gradient"))
bq <- summary(reg2Results1)[['bootstrapQuantile']]
omxCheckCloseEnough(cor(apply(bq,1,diff),
                        regResults1$output$standardErrors), 1.0, .2)
