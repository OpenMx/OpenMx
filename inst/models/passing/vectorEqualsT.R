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


# Purpose: Generate a script that causes an error by simply switching the fit function
#  from vector=FALSE to vector=TRUE.
#  The source of the script is Tim Bates' post here
# 	http://openmx.ssri.psu.edu/thread/861

require(OpenMx)

#Simulate Data
set.seed(200)
rs = .5
nSubs = 1000
selVars <- c('X','Y')
nVar = length(selVars)
xy <- mvtnorm::rmvnorm (nSubs, c(0,0), matrix(c(1,rs,rs,1),2,2))
testData <- data.frame(xy) 
names(testData) <- selVars
 
m1 <- mxModel("vector_is_FALSE", 
	mxMatrix(name = "expCov", type = "Symm", nrow = nVar, ncol = nVar, free = T, values = var(testData)),
	mxMatrix(name = "expMean", type = "Full", nrow = 1, ncol = nVar, free = T, values = 0),
	mxExpectationNormal(covariance = "expCov", means = "expMean", dimnames = selVars),
	mxFitFunctionML(vector = FALSE),
	mxData(observed = testData, type = "raw")
#  ,mxComputeGradientDescent(verbose=1L)
)
m1 <- mxRun(m1)

trueParam <- c(vech(cov(testData)), colMeans(testData))
omxCheckCloseEnough(omxGetParameters(m1), trueParam, 0.01)
 
# Now: Switch to vector
m2 <- mxModel(m1, mxFitFunctionML(vector = TRUE), name = "vector_is_TRUE")
m2 <- omxCheckWarning(mxRun(m2, suppressWarnings=TRUE), "vector_is_TRUE.fitfunction does not evaluate to a 1x1 matrix. Fixing model by adding mxAlgebra(-2*sum(log(vector_is_TRUE.fitfunction)), 'm2ll'), mxFitFunctionAlgebra('m2ll')")
omxCheckCloseEnough(omxGetParameters(m2), trueParam, 0.01)
