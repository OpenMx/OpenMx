#
#   Copyright 2007-2012 The OpenMx Project
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


# -----------------------------------------------------------------------------
# Program: BivariateSaturated_MatrixRawCholesky.R  
# Author: Hermine Maes
# Date: 2009.08.01 
#
# ModelType: Saturated
# DataType: Continuous
# Field: None
#
# Purpose: 
#      Bivariate Saturated model to estimate means and (co)variances 
#      using Cholesky Decomposition
#      Matrix style model input - Raw data input
#
# RevisionHistory:
#      Hermine Maes -- 2009.10.08 updated & reformatted
#      Ross Gore -- 2011.06.15 added Model, Data & Field metadata     
# -----------------------------------------------------------------------------

require(OpenMx)
require(MASS)
# Load Libraries
# -----------------------------------------------------------------------------

set.seed(200)
rs=.5
xy <- mvrnorm (1000, c(0,0), matrix(c(1,rs,rs,1),2,2))
testData <- xy
selVars <- c("X","Y")
dimnames(testData) <- list(NULL, selVars)
summary(testData)
cov(testData)
# Simulate Data
# -----------------------------------------------------------------------------

bivSatModel6 <- mxModel("bivSat6",
    mxMatrix(
		type="Full", 
		nrow=2, 
		ncol=2, 
		free=c(TRUE,TRUE,FALSE,TRUE), 
		values=c(1,.2,0,1), 
		name="Chol"
	),
	mxMatrix(
		type="Full", 
		nrow=1, 
		ncol=2, 
		free=TRUE, 
		values=c(0,0), 
		name="expMean"
	),
    mxAlgebra(
		expression=Chol %*% t(Chol), 
		name="expCov"
	),
	mxData(
		observed=testData, 
		type="raw"
	),
	mxFIMLObjective(
		covariance="expCov",
		means="expMean",
		dimnames=selVars
	)
)

bivSatFit6 <- mxRun(bivSatModel6)
EM <- mxEval(expMean, bivSatFit6)
EC <- mxEval(expCov, bivSatFit6)
LL <- mxEval(objective,bivSatFit6)
# example 6: Saturated Model with Raw Data and Matrix-Style Input
# -----------------------------------------------------------------------------

Mx.EM <- matrix(c(0.03211188, -0.004889211),1,2)
Mx.EC <- matrix(c(1.0092891, 0.4813504, 0.4813504, 0.9935366),2,2)
Mx.LL <- 5415.772
# Mx answers hard-coded
# -------------------------------------
# example Mx..2: Saturated Model with Raw Data
# -----------------------------------------------------------------------------

omxCheckCloseEnough(LL,Mx.LL,.001)
omxCheckCloseEnough(EC,Mx.EC,.001)
omxCheckCloseEnough(EM,Mx.EM,.001)
# 6: RawMat Cholesky
# -------------------------------------
# Compare OpenMx results to Mx results 
# (LL: likelihood; EC: expected covariance, EM: expected means)
# -----------------------------------------------------------------------------
