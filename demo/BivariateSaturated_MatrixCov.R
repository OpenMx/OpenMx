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
# Program: BivariateSaturated_MatrixCov.R  
# Author: Hermine Maes
# Date: 2009.08.01 
#
# ModelType: Saturated
# DataType: Continuous
# Field: None
#
# Purpose: 
#      Bivariate Saturated model to estimate means and (co)variances
#      Matrix style model input - Covariance matrix data input
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
covData <- cov(testData)
# Simulate Data
# -----------------------------------------------------------------------------


bivSatModel3 <- mxModel("bivSat3",
    mxMatrix(
        type="Symm", 
        nrow=2, 
        ncol=2, 
        free=T, 
        values=c(1,.5,1), 
        name="expCov"
    ),
    mxData(
        observed=covData, 
        type="cov", 
        numObs=1000 
    ),
    mxMLObjective(
        covariance="expCov",
        dimnames=selVars
    )
)

bivSatFit3 <- mxRun(bivSatModel3)
EC3 <- mxEval(expCov, bivSatFit3)
LL3 <- mxEval(objective,bivSatFit3)
bivSatSummary3 <- summary(bivSatFit3) 
SL3 <- bivSatSummary3$SaturatedLikelihood
omxCheckEquals(bivSatSummary3$observedStatistics, nrow(covData) * (ncol(covData) + 1) / 2)
Chi3 <- LL3-SL3
# example 3: Saturated Model with Cov Matrices and Matrix-Style Input
# -----------------------------------------------------------------------------


bivSatModel3m <- mxModel("bivSat3m",
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
        observed=cov(testData), 
        type="cov", 
        numObs=1000, 
        means=colMeans(testData) 
    ),
    mxMLObjective(
        covariance="expCov",
        means="expMean",
        dimnames=selVars
    )
)

bivSatFit3m <- mxRun(bivSatModel3m)
EM3m <- mxEval(expMean, bivSatFit3m)
EC3m <- mxEval(expCov, bivSatFit3m)
LL3m <- mxEval(objective,bivSatFit3m)
SL3m <- summary(bivSatFit3m)$SaturatedLikelihood
Chi3m <- LL3m-SL3m
# example 3m: Saturated Model with Cov Matrices & Means and Matrix-Style Input
# -----------------------------------------------------------------------------



Mx.EC1 <- matrix(c(1.0102951, 0.4818317, 0.4818317, 0.9945329),2,2)
Mx.LL1 <- -2.258885e-13
# example Mx..1: Saturated Model with 
# Cov Matrices
# -------------------------------------


Mx.EM1m <- matrix(c(0.03211648, -0.004883811),1,2)
Mx.EC1m <- matrix(c(1.0102951, 0.4818317, 0.4818317, 0.9945329),2,2)
Mx.LL1m <- -5.828112e-14
# example Mx..1m: Saturated Model with 
# Cov Matrices & Means
# -------------------------------------
# Mx answers hard-coded
# -----------------------------------------------------------------------------


omxCheckCloseEnough(Chi3,Mx.LL1,.001)
omxCheckCloseEnough(EC3,Mx.EC1,.001)
# 3:CovMat
# -------------------------------------

omxCheckCloseEnough(Chi3m,Mx.LL1m,.001)
omxCheckCloseEnough(EC3m,Mx.EC1m,.001)
omxCheckCloseEnough(EM3m,Mx.EM1m,.001)
# 3m:CovMPat 
# -------------------------------------
# Compare OpenMx results to Mx results 
# (LL: likelihood; EC: expected covariance, EM: expected means)
# -----------------------------------------------------------------------------
