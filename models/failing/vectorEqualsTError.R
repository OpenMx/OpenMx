#
#   Copyright 2007-2014 The OpenMx Project
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

#------------------------------------------------------------------------------
# Author: Michael D. Hunter
# Date: 2014.03.04
# Filename: vectorEqualsTError.R
# Purpose: Generate a script that causes an error by simply switching the fit function
#  from vector=FALSE to vector=TRUE.
#  The source of the script is Tim Bates' post here
# 	http://openmx.psyc.virginia.edu/thread/861



#------------------------------------------------------------------------------

require(OpenMx)
require(MASS)
#Simulate Data
set.seed(200)
rs = .5
nSubs = 1000
selVars <- c('X','Y')
nVar = length(selVars)
xy <- mvrnorm (nSubs, c(0,0), matrix(c(1,rs,rs,1),2,2))
testData <- data.frame(xy) 
names(testData) <- selVars
cov(testData)
 
m1 <- mxModel("vector_is_FALSE", 
	mxMatrix(name = "expCov", type = "Symm", nrow = nVar, ncol = nVar, free = T, values = var(testData)),
	mxMatrix(name = "expMean", type = "Full", nrow = 1, ncol = nVar, free = T, values = 0),
	mxExpectationNormal(covariance = "expCov", means = "expMean", dimnames = selVars),
	mxFitFunctionML(vector = FALSE),
	mxData(observed = testData, type = "raw")
)
m1 <- mxRun(m1)

trueParam <- c(vech(cov(testData)), colMeans(testData))
omxCheckCloseEnough(omxGetParameters(m1), trueParam, 0.01)
 
# Now: Switch to vector
m1 <- mxModel(m1, mxFitFunctionML(vector = TRUE), name = "vector_is_TRUE")
omxCheckError(mxRun(m1), 
	paste("The top level fitfunction in model", m1$name, "does not",
	"evaluate to a 1x1 matrix.  Consider aggregating by adding the following to your model:",
	"mxAlgebra(-2*sum(log(YourModelName.fitfunction))), 'm2ll'), mxFitFunctionAlgebra('m2ll')"))

# what we get
#                  name  Estimate Std.Error
# 1  vector.expCov[1,1]    0.4705   9.0e+07
# 2  vector.expCov[1,2]    0.6148   1.4e+08
# 3  vector.expCov[2,2]    1.0314   2.1e+08
# what it should be...
cov(testData)
# XX 0.9945328
# XY 0.4818317 
# YY 1.0102951
