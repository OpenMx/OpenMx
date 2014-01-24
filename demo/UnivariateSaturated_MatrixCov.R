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
# Program: UnivariateSaturated_MatrixCov.R  
# Author: Hermine Maes
# Date: 2009.08.01 
#
# ModelType: Saturated
# DataType: Simulated
# Field: None
#
# Purpose: 
#      Univariate Saturated model to estimate means and variances
#      Matrix style model input - Covariance matrix data input
#
# RevisionHistory:
#      Hermine Maes -- 2009.10.08 updated & reformatted
#      Ross Gore -- 2011.06.06 added Model, Data & Field metadata
# -----------------------------------------------------------------------------

require(OpenMx)
# Load Library
# -----------------------------------------------------------------------------

set.seed(100)
x <- rnorm (1000, 0, 1)
testData <- as.matrix(x)
selVars <- c("X")
dimnames(testData) <- list(NULL, selVars)
summary(testData)
colMeans(testData)
var(testData)
# Simulate Data
# -----------------------------------------------------------------------------


univSatModel3 <- mxModel("univSat3",
    mxMatrix(
        type="Symm", 
        nrow=1, 
        ncol=1, 
        free=T, 
        values=1, 
        name="expCov"
    ),
    mxData(
        observed=var(testData), 
        type="cov", 
        numObs=1000 
    ),
    mxMLObjective(
        covariance="expCov",
        dimnames=selVars
    )
)

univSatFit3 <- mxRun(univSatModel3)
EC3 <- mxEval(expCov, univSatFit3)
LL3 <- mxEval(objective, univSatFit3)
SL3 <- univSatFit3@output$SaturatedLikelihood
Chi3 <- LL3-SL3
# example 3: Saturated Model with Cov Matrices and Matrix-Style Input
# -----------------------------------------------------------------------------

univSatModel3m <- mxModel("univSat3m",
    mxMatrix(
        type="Symm", 
        nrow=1, 
        ncol=1, 
        free=T, 
        values=1, 
        name="expCov"
    ),
    mxMatrix(
        type="Full", 
        nrow=1, 
        ncol=1, 
        free=T, 
        values=0, 
        name="expMean"
    ),
    mxData(
        observed=var(testData), 
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

univSatFit3m <- mxRun(univSatModel3m)
EM3m <- mxEval(expMean, univSatFit3m)
EC3m <- mxEval(expCov, univSatFit3m)
LL3m <- mxEval(objective, univSatFit3m);
SL3m <- univSatFit3m@output$SaturatedLikelihood
Chi3m <- LL3m-SL3m
# example 3m: Saturated Model with Cov Matrices & Means 
# and Matrix-Style Input
# -----------------------------------------------------------------------------



Mx.EC1 <-  1.062112
Mx.LL1 <- -1.474434e-17
# example Mx..1: Saturated Model with 
# Cov Matrices
# -------------------------------------

Mx.EM1m <- 0.01680509
Mx.EC1m <- 1.062112
Mx.LL1m <- -1.108815e-13
# example Mx..1m: Saturated Model with 
# Cov Matrices & Means
# -------------------------------------
# Mx answers hard-coded
# -----------------------------------------------------------------------------


##3:CovMat
# -------------------------------------
omxCheckCloseEnough(Chi3,Mx.LL1,.001)
omxCheckCloseEnough(EC3,Mx.EC1,.001)
#3m:CovMPat 
# -------------------------------------
omxCheckCloseEnough(Chi3m,Mx.LL1m,.001)
omxCheckCloseEnough(EC3m,Mx.EC1m,.001)
omxCheckCloseEnough(EM3m,Mx.EM1m,.001)
# Compare OpenMx results to Mx results 
# (LL: likelihood; EC: expected covariance, EM: expected means)
# -----------------------------------------------------------------------------
