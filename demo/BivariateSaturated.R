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
# Program: BivariateSaturated.R  
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
#      Two matrix styles - Two data styles
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
xy <- mvtnorm::rmvnorm (1000, c(0,0), matrix(c(1,rs,rs,1),2,2))
testData <- xy
testData <- testData[, order(apply(testData, 2, var))[2:1]] #put the data columns in order from largest to smallest variance
# Note: Users do NOT have to re-order their data columns.  This is only to make data generation the same on different operating systems: to fix an inconsistency with the mvtnorm::rmvnorm function in the MASS package.
selVars <- c('X','Y')
dimnames(testData) <- list(NULL, selVars)
summary(testData)
cov(testData)
# Simulate Data
# -----------------------------------------------------------------------------

bivSatModel1 <- mxModel("bivSat1",
    manifestVars= selVars,
    mxPath(
        from=c("X", "Y"), 
        arrows=2, 
        free=T, 
        values=1, 
        lbound=.01, 
        labels=c("varX","varY")
    ),
    mxPath(
        from="X", 
        to="Y", 
        arrows=2, 
        free=T, 
        values=.2, 
        lbound=.01, 
        labels="covXY"
    ),
    mxData(
        observed=cov(testData), 
        type="cov", 
        numObs=1000 
    ),
    type="RAM"
    )
bivSatFit1 <- mxRun(bivSatModel1)
EC1 <- mxEval(S, bivSatFit1)
LL1 <- mxEval(fitfunction, bivSatFit1)
SL1 <- summary(bivSatFit1)$SaturatedLikelihood
Chi1 <- LL1-SL1
# example 1: Saturated Model with Cov Matrices and Path-Style Input
# -----------------------------------------------------------------------------

bivSatModel1m <- mxModel("bivSat1m",
    manifestVars= selVars,
    mxPath(
        from=c("X", "Y"), 
        arrows=2, 
        free=T, 
        values=1, 
        lbound=.01, 
        labels=c("varX","varY")
    ),
    mxPath(
        from="X", 
        to="Y", 
        arrows=2, 
        free=T, 
        values=.2, 
        lbound=.01, 
        labels="covXY"
    ),
    mxPath(
        from="one", 
        to=c("X", "Y"), 
        arrows=1, 
        free=T, 
        values=.01, 
        labels=c("meanX","meanY")
    ),
    mxData(
        observed=cov(testData), 
        type="cov", 
        numObs=1000, 
        means=colMeans(testData)
    ),
    type="RAM"
    )
bivSatFit1m <- mxRun(bivSatModel1m)
EM1m <- mxEval(M, bivSatFit1m)
EC1m <- mxEval(S, bivSatFit1m)
LL1m <- mxEval(fitfunction, bivSatFit1m)
SL1m <- summary(bivSatFit1m)$SaturatedLikelihood
Chi1m <- LL1m-SL1m
# example 1m: Saturated Model with Cov Matrices & Means and Path-Style Input
# -----------------------------------------------------------------------------

bivSatModel2 <- mxModel("bivSat2",
    manifestVars= selVars,
    mxPath(
        from=c("X", "Y"), 
        arrows=2, 
        free=T, 
        values=1, 
        lbound=.01, 
        labels=c("varX","varY")
    ),
    mxPath(
        from="X", 
        to="Y", 
        arrows=2, 
        free=T, 
        values=.2, 
        lbound=.01, 
        labels="covXY"
    ),
    mxPath(from="one", 
    	to=c("X", "Y"),
        arrows=1, 
        free=T), 
    mxData(
        observed=testData, 
        type="raw", 
    ),
    type="RAM"
    )
bivSatFit2 <- mxRun(bivSatModel2)
EM2 <- mxEval(M, bivSatFit2)
EC2 <- mxEval(S, bivSatFit2)
LL2 <- mxEval(fitfunction, bivSatFit2)
SL2 <- summary(bivSatFit1)$SaturatedLikelihood
Chi2 <- LL2-SL2
# example 2: Saturated Model with Raw Data and Path input
# -----------------------------------------------------------------------------

bivSatModel2s <- mxModel(bivSatModel1,
    mxData(
        observed=testData, 
        type="raw", 
    ),
    mxPath(from="one", 
    	to=c("X", "Y"),
        arrows=1, 
        free=T),     
    name = "bivSat2s"
    )
bivSatFit2s <- mxRun(bivSatModel2s)
EM2s <- mxEval(M, bivSatFit2s)
EC2s <- mxEval(S, bivSatFit2s)
LL2s <- mxEval(fitfunction, bivSatFit2s)
# example 2s: Saturated Model with Raw Data and Path input built upon Cov/Means version
# -----------------------------------------------------------------------------

bivSatModel3 <- mxModel("bivSat3",
    mxMatrix(
        type="Symm", 
        nrow=2, 
        ncol=2, 
        free=T, 
        values=c(1,.5,1), 
        dimnames=list(selVars,selVars), 
        name="expCov"
    ),
    mxData(
        observed=cov(testData), 
        type="cov", 
        numObs=1000 
    ),
    mxExpectationNormal(
        covariance="expCov"
    ),
    mxFitFunctionML()
    )
bivSatFit3 <- mxRun(bivSatModel3)
EC3 <- mxEval(expCov, bivSatFit3)
LL3 <- mxEval(fitfunction,bivSatFit3)
SL3 <- summary(bivSatFit3)$SaturatedLikelihood
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
    mxExpectationNormal(
        covariance="expCov",
        means="expMean",
        dimnames=selVars
    ),
    mxFitFunctionML()
    )
bivSatFit3m <- mxRun(bivSatModel3m)
EM3m <- mxEval(expMean, bivSatFit3m)
EC3m <- mxEval(expCov, bivSatFit3m)
LL3m <- mxEval(fitfunction,bivSatFit3m)
SL3m <- summary(bivSatFit3m)$SaturatedLikelihood
Chi3m <- LL3m-SL3m
# example 3m: Saturated Model with Cov Matrices & Means and Matrix-Style Input
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
mxExpectationNormal(
    covariance="expCov",
    means="expMean",
    dimnames=selVars
),
mxFitFunctionML()
)
bivSatFit4 <- mxRun(bivSatModel4)
EM4 <- mxEval(expMean, bivSatFit4)
EC4 <- mxEval(expCov, bivSatFit4)
LL4 <- mxEval(fitfunction,bivSatFit4)
# examples 4: Saturated Model with Raw Data and Matrix-Style Input
# -----------------------------------------------------------------------------

bivSatModel5 <- mxModel("bivSat5",
    mxMatrix(
        type="Full", 
        nrow=2, 
        ncol=2, 
        free=c(T,T,F,T), 
        values=c(1,.2,0,1), 
        name="Chol"
    ),
    mxAlgebra(
        Chol %*% t(Chol), 
        name="expCov", 
        dimnames=list(selVars,selVars)
    ),
    mxMatrix(
        type="Full", 
        nrow=1, 
        ncol=2, 
        free=T, 
        values=c(0,0), 
        dimnames=list(NULL,selVars), 
        name="expMean"
    ),
    mxData(
        observed=cov(testData), 
        type="cov", 
        numObs=1000 
    ),
    mxExpectationNormal(
        covariance="expCov"
    ),
    mxFitFunctionML()
    )
bivSatFit5 <- mxRun(bivSatModel5)
EC5 <- mxEval(expCov, bivSatFit5)
LL5 <- mxEval(fitfunction,bivSatFit5)
SL5 <- summary(bivSatFit5)$SaturatedLikelihood
Chi5 <- LL5-SL5
# example 5: Saturated Model with Cov Matrices and Matrix-Style Input
# -----------------------------------------------------------------------------

bivSatModel5m <- mxModel("bivSat5m",
    mxMatrix(
        type="Full", 
        nrow=2, 
        ncol=2, 
        free=c(T,T,F,T), 
        values=c(1,.2,0,1), 
        name="Chol"
    ),
    mxAlgebra(
        Chol %*% t(Chol), 
        name="expCov", 
        dimnames=list(selVars,selVars)
    ),
    mxMatrix(
        type="Full", 
        nrow=1, 
        ncol=2, 
        free=T, 
        values=c(0,0), 
        dimnames=list(NULL,selVars), 
        name="expMean"
    ),
    mxData(
        observed=cov(testData), 
        type="cov", 
        numObs=1000, 
        means=colMeans(testData) 
    ),
    mxExpectationNormal(
        covariance="expCov",
        means="expMean"
    ),
    mxFitFunctionML()
    )
bivSatFit5m <- mxRun(bivSatModel5m)
EM5m <- mxEval(expMean, bivSatFit5m)
EC5m <- mxEval(expCov, bivSatFit5m)
LL5m <- mxEval(fitfunction,bivSatFit5m);
SL5m <- summary(bivSatFit5m)$SaturatedLikelihood
Chi5m <- LL5m-SL5m
# example 5m: Saturated Model with Cov Matrices & Means and Matrix-Style Input
# -----------------------------------------------------------------------------

bivSatModel6 <- mxModel("bivSat6",
    mxMatrix(
        type="Full", 
        nrow=2, 
        ncol=2, 
        free=c(T,T,F,T), 
        values=c(1,.2,0,1), 
        name="Chol"
    ),
    mxAlgebra(
        Chol %*% t(Chol), 
        name="expCov", 
        dimnames=list(selVars,selVars)
    ),
    mxMatrix(
        type="Full", 
        nrow=1, 
        ncol=2, 
        free=T, 
        values=c(0,0), 
        dimnames=list(NULL,selVars), 
        name="expMean"
    ),
    mxData(
        observed=testData, 
        type="raw", 
    ),
    mxExpectationNormal(
        covariance="expCov",
        means="expMean"
    ),
    mxFitFunctionML()
    )
bivSatFit6 <- mxRun(bivSatModel6)
EM6 <- mxEval(expMean, bivSatFit6)
EC6 <- mxEval(expCov, bivSatFit6)
LL6 <- mxEval(fitfunction,bivSatFit6)
# examples 6: Saturated Model with RawData and Matrix-Style Input
# -----------------------------------------------------------------------------


omxCheckCloseEnough(Chi1, 0,.01)
omxCheckCloseEnough(c(EC1), c(1.0656, 0.4751, 0.4751, 0.9291),.001)
#1:CovPat
# -------------------------------------

omxCheckCloseEnough(Chi1m, 0,.01)
omxCheckCloseEnough(c(EC1m), c(1.0656, 0.4751, 0.4751, 0.9291),.001)
omxCheckCloseEnough(c(EM1m), c(0.0582, 0.0063), .001)
#1m:CovMPat 
# -------------------------------------

omxCheckCloseEnough(LL2, 5407.037,.01)
omxCheckCloseEnough(c(EC2), c(1.0656, 0.4752, 0.4752, 0.9292),.001)
omxCheckCloseEnough(c(EM2), c(0.0582, 0.0063), .001)
#2:RawPat 
# -------------------------------------

omxCheckCloseEnough(LL2s, 5407.037, .01)
omxCheckCloseEnough(c(EC2s), c(1.0656, 0.4752, 0.4752, 0.9292),.001)
omxCheckCloseEnough(c(EM2s), c(0.0582, 0.0063),.001)
#2:RawSPat
# -------------------------------------

omxCheckCloseEnough(Chi3, 0,.01)
omxCheckCloseEnough(c(EC3),c(1.0656, 0.4752, 0.4752, 0.9292),.001)
#3:CovMat
# -------------------------------------

omxCheckCloseEnough(Chi3m,0,.01)
omxCheckCloseEnough(c(EC3m),c(1.0656, 0.4752, 0.4752, 0.9292),.001)
omxCheckCloseEnough(c(EM3m),c(0.058, 0.006),.001)
#3m:CovMPat 
# -------------------------------------

omxCheckCloseEnough(LL4,5407.037,.01)
omxCheckCloseEnough(c(EC4),c(1.066, 0.475, 0.475, 0.929),.001)
omxCheckCloseEnough(c(EM4),c(0.058, 0.006),.001)
#4:RawMat
# -------------------------------------

omxCheckCloseEnough(Chi5,0,.01)
omxCheckCloseEnough(c(EC5),c(1.0656, 0.4752, 0.4752, 0.9292),.001)
#5:CovMat Cholesky
# -------------------------------------

omxCheckCloseEnough(Chi5m,0,.01)
omxCheckCloseEnough(c(EC5m),c(1.0656, 0.4752, 0.4752, 0.9292),.001)
omxCheckCloseEnough(c(EM5m),c(0.058, 0.006),.001)
#5m:CovMPat Cholesky
# -------------------------------------


omxCheckCloseEnough(LL6,5407.037,.01)
omxCheckCloseEnough(c(EC6),c(1.066, 0.475, 0.475, 0.929),.001)
omxCheckCloseEnough(c(EM6),c(0.058, 0.006),.001)
#6:RawMat Cholesky
# -------------------------------------
# Compare OpenMx results to Mx results 
# (LL: likelihood; EC: expected covariance, EM: expected means)
# -----------------------------------------------------------------------

