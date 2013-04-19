#
#   Copyright 2007-2013 The OpenMx Project
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
# -----------------------------------------------------------------------------

require(OpenMx)
require(MASS)
# Load Library
# -----------------------------------------------------------------------------

set.seed(200)
rs=.5
xy <- mvrnorm (1000, c(0,0), matrix(c(1,rs,rs,1),2,2))
testData <- xy
testData <- testData[, order(apply(testData, 2, var))[2:1]] #put the data columns in order from largest to smallest variance
selVars <- c('X','Y')
dimnames(testData) <- list(NULL, selVars)
summary(testData)
cov(testData)
# Simulate Data: two standardized variables X & Y with correlation of .5
# -----------------------------------------------------------------------------


bivCorModel <- mxModel("bivCor",
    mxMatrix(
        type="Full", 
        nrow=1, 
        ncol=2, 
        free=TRUE, 
        values=c(0,0), 
        name="expMean"
    ), 
    mxMatrix(
        type="Lower", 
        nrow=2, 
        ncol=2, 
        free=TRUE,
        values=.5, 
        name="Chol"
    ), 
    mxAlgebra(
        expression=Chol %*% t(Chol), 
        name="expCov", 
    ), 
    mxData(
        observed=testData, 
        type="raw"
    ), 
    mxExpectationNormal(
        covariance="expCov", 
        means="expMean",
        dimnames=selVars),
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
    mxMatrix(
        type="Diag", 
        nrow=2, 
        ncol=2, 
        free=TRUE, 
        name="Chol"
    )
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


Mx.EM <- matrix(c(0.03211656, -0.004883885), 1, 2)
Mx.EC <- matrix(c(1.0092853, 0.4813504, 0.4813504, 0.9935390), 2, 2)
Mx.LL <- 5415.772
# Mx Answers of Saturated Model Hard-coded
# -----------------------------------------------------------------------------

omxCheckCloseEnough(LL,Mx.LL,.001)
omxCheckCloseEnough(EC,Mx.EC,.001)
omxCheckCloseEnough(EM,Mx.EM,.001)
# Compare OpenMx Results to Mx Results 
# LL: likelihood; EC: expected covariance, EM: expected means
# -----------------------------------------------------------------------------
