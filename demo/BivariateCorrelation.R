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
# Program: BivariateCorrelation.R  
# Author: Hermine Maes
# Date: 2009.08.01 
#
# ModelType: Saturated
# DataType: Continuous
# Field: None
#
# Purpose: 
#      Optimization Example in OpenMx: Testing significance of correlation
#
# RevisionHistory:
#      Hermine Maes -- 2009.10.08 updated & reformatted
#      Ross Gore -- 2011.06.15 added Model, Data & Field metadata
#      Mike Hunter -- 2013.09.16 nudged starting values of second model varainces away from zero
#      Hermine Maes -- 2014.11.02 piecewise specification
# -----------------------------------------------------------------------------

require(OpenMx)
require(MASS)
# Load Library
# -----------------------------------------------------------------------------

set.seed(200)
rs=.5
xy <- mvtnorm::rmvnorm (1000, c(0,0), matrix(c(1,rs,rs,1),2,2))
testData <- xy
testData <- testData[, order(apply(testData, 2, var))[2:1]] #put the data columns in order from largest to smallest variance
# Note: Users do NOT have to re-order their data columns.  This is only to make data generation the same on different operating systems: to fix an inconsistency with the mvtnorm::rmvnorm function in the MASS package.
selVars <- c('X','Y')
dimnames(testData) <- list(NULL, selVars)
summary(testData)
cov(testData)
# Simulate Data: two standardized variables X & Y with correlation of .5
# -----------------------------------------------------------------------------


bivCorModel <- mxModel("bivCor",
    mxMatrix( type="Full", nrow=1, ncol=2, free=TRUE, values=c(0,0), name="expMean" ), 
    mxMatrix( type="Lower", nrow=2, ncol=2, free=TRUE, values=.5, name="Chol" ), 
    mxAlgebra( expression=Chol %*% t(Chol), name="expCov"), 
    mxData( observed=testData, type="raw" ), 
    mxExpectationNormal( covariance="expCov", means="expMean", dimnames=selVars),
    mxFitFunctionML()
    )
# Fit Saturated Model with Raw Data and Matrix-style Input
# -----------------------------------------------------------------------------


bivCorFit <- mxRun(bivCorModel)
EM <- mxEval(expMean, bivCorFit)
EC <- mxEval(expCov, bivCorFit)
LL <- mxEval(fitfunction, bivCorFit)
# Run Model and Generate Output
# -----------------------------------------------------------------------------


bivCorModelSub <-mxModel(bivCorModel,
    mxMatrix( type="Diag", nrow=2, ncol=2, free=TRUE,
        values=.2, # Note: to test optimizer for robustness to bad starting values, change to 0.
        name="Chol" )
)
# Specify SubModel testing Covariance=Zero
# -----------------------------------------------------------------------------


bivCorFitSub <- mxRun(bivCorModelSub)
EMs <- mxEval(expMean, bivCorFitSub)
ECs <- mxEval(expCov, bivCorFitSub)
LLs <- mxEval(fitfunction, bivCorFitSub)
Chi= LLs-LL;
LRT= rbind(LL,LLs,Chi); LRT
# Run Model and Generate Output
# -----------------------------------------------------------------------------


omxCheckCloseEnough(LL, 5407.036, .001)
omxCheckCloseEnough(c(EC), c(1.0656, 0.4752, 0.4752, 0.9292), .001)
omxCheckCloseEnough(c(EM), c(0.058, 0.006), .001)
# Compare OpenMx Results to Mx Results 
# LL: likelihood; EC: expected covariance, EM: expected means
# -----------------------------------------------------------------------------
