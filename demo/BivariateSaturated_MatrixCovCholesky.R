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


# -----------------------------------------------------------------------
# Program: BivariateSaturated_MatrixCovCholesky.R  
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
#      Matrix style model input - Covariance matrix data input
#
# RevisionHistory:
#      Hermine Maes -- 2009.10.08 updated & reformatted
#      Ross Gore -- 2011.06.15 added Model, Data & Field metadata
# -----------------------------------------------------------------------

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


bivSatModel5 <- mxModel("bivSat5",
    mxMatrix(
        type="Lower", 
        nrow=2, 
        ncol=2, 
        free=T, 
        values=.5, 
        name="Chol"
    ),
    mxAlgebra(
        expression=Chol %*% t(Chol), 
        name="expCov" 
    ),
    mxData(
        observed=cov(testData), 
        type="cov", 
        numObs=1000 
    ),
    mxMLObjective(
        covariance="expCov",
        dimnames=selVars
    )
)

bivSatFit5 <- mxRun(bivSatModel5)
EC5 <- mxEval(expCov, bivSatFit5)
LL5 <- mxEval(objective,bivSatFit5)
SL5 <- summary(bivSatFit5)$SaturatedLikelihood
Chi5 <- LL5-SL5
# example 5: Saturated Model with Cov Matrices and Matrix-Style Input
# -----------------------------------------------------------------------------


bivSatModel5m <- mxModel("bivSat5m",
    mxMatrix(
        type="Lower", 
        nrow=2, 
        ncol=2, 
        free=T, 
        values=.5, 
        name="Chol"
    ),
    mxAlgebra(
        expression=Chol %*% t(Chol), 
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

bivSatFit5m <- mxRun(bivSatModel5m)
EM5m <- mxEval(expMean, bivSatFit5m)
EC5m <- mxEval(expCov, bivSatFit5m)
LL5m <- mxEval(objective,bivSatFit5m);
SL5m <- summary(bivSatFit5m)$SaturatedLikelihood
Chi5m <- LL5m-SL5m
# example 5m: Saturated Model with Cov Matrices & Means and Matrix-Style Input
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



omxCheckCloseEnough(Chi5,Mx.LL1,.001)
omxCheckCloseEnough(EC5,Mx.EC1,.001)
# 5:CovMat Cholesky
# -------------------------------------

omxCheckCloseEnough(Chi5m,Mx.LL1m,.001)
omxCheckCloseEnough(EC5m,Mx.EC1m,.001)
omxCheckCloseEnough(EM5m,Mx.EM1m,.001)
# 5m:CovMPat Cholesky
# -------------------------------------
# Compare OpenMx results to Mx results 
# (LL: likelihood; EC: expected covariance, EM: expected means)
# -----------------------------------------------------------------------------
