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


# -----------------------------------------------------------------------------
# Program: BivariateSaturated_MatrixRaw.R  
# Author: Hermine Maes
# Date: 2009.08.01 
#
# ModelType: Saturated
# DataType: Continuous
# Field: None
#
# Purpose:
#      Bivariate Saturated model to estimate means and (co)variances
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
xy <- mvtnorm::rmvnorm (1000, c(0,0), matrix(c(1,rs,rs,1),2,2))
testData <- xy
testData <- testData[, order(apply(testData, 2, var))[2:1]] #put the data columns in order from largest to smallest variance
# Note: Users do NOT have to re-order their data columns.  This is only to make data generation the same on different operating systems: to fix an inconsistency with the mvtnorm::rmvnorm function in the MASS package.
selVars <- c("X","Y")
dimnames(testData) <- list(NULL, selVars)
summary(testData)
cov(testData)
# Simulate Data
# -----------------------------------------------------------------------------

bivSatModel4 <- mxModel("bivSat4",
	mxMatrix(
	    type="Symm", 
	    nrow=2, 
	    ncol=2, 
	    free=T, 
	    values=c(1,.5,1), 
	    name="expCov"
	),
	mxMatrix(
	    type="Full", 
	    nrow=1, 
	    ncol=2, 
	    free=T, 
	    values=c(0,0), 
	    name="expMean"
	),
	mxData(
	    observed=testData, 
	    type="raw", 
	),
	mxFitFunctionML(),mxExpectationNormal(
	    covariance="expCov",
	    means="expMean",
	    dimnames=selVars
	)
)

bivSatFit4 <- mxRun(bivSatModel4)
bivSatSummary4 <- summary(bivSatFit4)
EM4 <- mxEval(expMean, bivSatFit4)
EC4 <- mxEval(expCov, bivSatFit4)
LL4 <- mxEval(objective,bivSatFit4)
omxCheckEquals(bivSatSummary4$observedStatistics, sum(!is.na(testData)))
# examples 4: Saturated Model with Raw Data and Matrix-Style Input
# -----------------------------------------------------------------------------

omxCheckCloseEnough(LL4,5407.037,.001)
omxCheckCloseEnough(c(EC4),c(1.066, 0.475, 0.475, 0.929),.001)
omxCheckCloseEnough(c(EM4),c(0.058, 0.006),.001)
# 4:RawMat
# -------------------------------------
# Compare OpenMx results to Mx results 
# (LL: likelihood; EC: expected covariance, EM: expected means)
# -----------------------------------------------------------------------------
