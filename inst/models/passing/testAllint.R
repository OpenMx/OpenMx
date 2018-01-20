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


# -----------------------------------------------------------------------
# Program: testAllint.R  
#  Author: Mike Neale and Tim Brick
#    Date: 2010-02-20
#
#  Simple tests for omxAllint().
# 
# -----------------------------------------------------------------------

require(OpenMx)

mxOption(NULL, 'mvnMaxPointsC', 100)

# Simple case: 1 var, 3 levels

nv <- 1
maxnt <- 3
nvminus1 <- nv - 1

testAllint1 <- mxModel("TestAllint",
		mxMatrix(type="Stand", nrow=nv, free=F, name="Cov"),
		mxMatrix(type="Full", nrow=1, ncol=nv, free=F, name="Means"),
		mxMatrix(type="Full", nrow=maxnt, ncol=nv, free=F, values = c(-Inf, 2.326, Inf), name="Thresh"),
		mxAlgebra(omxAllInt(Cov, Means, Thresh), name="testAllint"))

testAllintFit1 <- mxRun(testAllint1)

#Test 2x3

nv <- 2
maxnt <- 3
nvminus1 <- nv - 1

testAllint2 <- mxModel("TestAllint",
		mxMatrix(type="Stand", nrow=nv, free=F, name="Cov"),
		mxMatrix(type="Full", nrow=1, ncol=nv, free=F, name="Means"),
		mxMatrix(type="Full", free=F, values = cbind(c(-Inf, 0, Inf), c(-Inf, 1.96, Inf)), name="Thresh"),
		mxAlgebra(omxAllInt(Cov, Means, Thresh), name="testAllint"))

testAllintFit2 <- mxRun(testAllint2)

#Test 2x4

nv <- 2
maxnt <- 4
nvminus1 <- nv - 1

testAllint3 <- mxModel("TestAllint",
		mxMatrix(type="Stand", nrow=nv, free=F, name="Cov"),
		mxMatrix(type="Full", nrow=1, ncol=nv, free=F, name="Means"),
		mxMatrix(type="Full", free=F, values = cbind(c(-Inf, 0, 1, Inf), c(-Inf, 1.96, 2.326, Inf)), name="Thresh"),
		mxAlgebra(omxAllInt(Cov, Means, Thresh), name="testAllint"))

testAllintFit3 <- mxRun(testAllint3)

#Test 2x4,1x5 with NAs

nv <- 3
maxnt <- 5
nvminus1 <- nv - 1

testAllint4 <- mxModel("TestAllint",
		mxMatrix(type="Stand", nrow=nv, free=F, name="Cov"),
		mxMatrix(type="Full", nrow=nv, ncol=1, free=F, name="Means"),
		mxMatrix(type="Full", free=F, values = cbind(c(-Inf, 0, 1, Inf, NA), c(-Inf, 1.96, 2.326, Inf, NA), c(-Inf, -1, 0, 1, Inf)), name="Thresh"),
		mxAlgebra(omxAllInt(Cov, Means, Thresh), name="testAllint"))

testAllintFit4 <- mxRun(testAllint4)

# Test different sizes of matrix
nv <- 3
maxnt <- 5
nvminus1 <- nv - 1

testAllint5 <- mxModel("TestAllint",
		mxMatrix(type="Stand", nrow=nv, free=F, name="Cov"),
		mxMatrix(type="Full", nrow=1, ncol=nv, free=F, name="Means"),
		mxMatrix(type="Full", free=F, values = cbind(c(-Inf, 0, 1, Inf), c(-Inf, 1.96, 2.326, Inf)), name="Thresh1"),
		mxMatrix(type="Full", free=F, values = cbind(c(-Inf, -1, 0, 1, Inf)), name="Thresh2"), 
		mxAlgebra(omxAllInt(Cov, Means, Thresh1, Thresh2), name="testAllint"))

testAllintFit5 <- mxRun(testAllint5)

#Test against Mx1 solutions
omxCheckCloseEnough(testAllintFit1[['testAllint']]$result, as.matrix(c(.99, .01)), 0.001)
omxCheckCloseEnough(testAllintFit2[['testAllint']]$result, as.matrix(c(0.4875,  0.0125, 0.4875, 0.0125)), 0.001)
omxCheckCloseEnough(testAllintFit3[['testAllint']]$result, as.matrix(c(0.4875,  0.0075, 0.0050, 0.3328, 0.0051, 0.0034, 0.1547, 0.0024, 0.0016)), 0.001)
omxCheckCloseEnough(testAllintFit4[['testAllint']]$result, as.matrix(c(7.7345E-02, 1.6641E-01, 1.6641E-01, 7.7345E-02, 1.1890E-03, 2.5581E-03, 2.5581E-03, 1.1890E-03, 7.9401E-04, 1.7083E-03, 1.7083E-03, 7.9401E-04, 5.2802E-02, 1.1360E-01, 1.1360E-01, 5.2802E-02, 8.1173E-04,  1.7464E-03, 1.7464E-03, 8.1173E-04, 5.4206E-04, 1.1662E-03, 1.1662E-03, 5.4206E-04, 2.4542E-02, 5.2802E-02, 5.2802E-02, 2.4542E-02, 3.7729E-04,  8.1173E-04, 8.1173E-04, 3.7729E-04, 2.5195E-04, 5.4206E-04, 5.4206E-04, 2.5195E-04)), 0.001)
omxCheckCloseEnough(testAllintFit5[['testAllint']]$result, as.matrix(c(7.7345E-02, 1.6641E-01, 1.6641E-01, 7.7345E-02, 1.1890E-03, 2.5581E-03, 2.5581E-03, 1.1890E-03, 7.9401E-04, 1.7083E-03, 1.7083E-03, 7.9401E-04, 5.2802E-02, 1.1360E-01, 1.1360E-01, 5.2802E-02, 8.1173E-04,  1.7464E-03, 1.7464E-03, 8.1173E-04, 5.4206E-04, 1.1662E-03, 1.1662E-03, 5.4206E-04, 2.4542E-02, 5.2802E-02, 5.2802E-02, 2.4542E-02, 3.7729E-04,  8.1173E-04, 8.1173E-04, 3.7729E-04, 2.5195E-04, 5.4206E-04, 5.4206E-04, 2.5195E-04)), 0.001)

omxCheckCloseEnough(mxEval(omxAllInt(Cov, Means, Thresh), testAllint1), as.matrix(c(.99, .01)), 0.001)
omxCheckCloseEnough(mxEval(omxAllInt(Cov, Means, Thresh), testAllint2), as.matrix(c(0.4875,  0.0125, 0.4875, 0.0125)), 0.001)
omxCheckCloseEnough(mxEval(omxAllInt(Cov, Means, Thresh), testAllint3), as.matrix(c(0.4875,  0.0075, 0.0050, 0.3328, 0.0051, 0.0034, 0.1547, 0.0024, 0.0016)), 0.001)
omxCheckCloseEnough(mxEval(omxAllInt(Cov, Means, Thresh), testAllint4), as.matrix(c(7.7345E-02, 1.6641E-01, 1.6641E-01, 7.7345E-02, 1.1890E-03, 2.5581E-03, 2.5581E-03, 1.1890E-03, 7.9401E-04, 1.7083E-03, 1.7083E-03, 7.9401E-04, 5.2802E-02, 1.1360E-01, 1.1360E-01, 5.2802E-02, 8.1173E-04,  1.7464E-03, 1.7464E-03, 8.1173E-04, 5.4206E-04, 1.1662E-03, 1.1662E-03, 5.4206E-04, 2.4542E-02, 5.2802E-02, 5.2802E-02, 2.4542E-02, 3.7729E-04,  8.1173E-04, 8.1173E-04, 3.7729E-04, 2.5195E-04, 5.4206E-04, 5.4206E-04, 2.5195E-04)), 0.001)
omxCheckCloseEnough(mxEval(omxAllInt(Cov, Means, Thresh1, Thresh2), testAllint5), as.matrix(c(7.7345E-02, 1.6641E-01, 1.6641E-01, 7.7345E-02, 1.1890E-03, 2.5581E-03, 2.5581E-03, 1.1890E-03, 7.9401E-04, 1.7083E-03, 1.7083E-03, 7.9401E-04, 5.2802E-02, 1.1360E-01, 1.1360E-01, 5.2802E-02, 8.1173E-04,  1.7464E-03, 1.7464E-03, 8.1173E-04, 5.4206E-04, 1.1662E-03, 1.1662E-03, 5.4206E-04, 2.4542E-02, 5.2802E-02, 5.2802E-02, 2.4542E-02, 3.7729E-04,  8.1173E-04, 8.1173E-04, 3.7729E-04, 2.5195E-04, 5.4206E-04, 5.4206E-04, 2.5195E-04)), 0.001)
