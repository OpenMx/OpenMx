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
# Program: UnivariateSaturated.R  
# Author: Hermine Maes
# Date: 2009.08.01 
#
# ModelType: Saturated
# DataType: Simulated
# Field: None
#
# Purpose: 
#      Univariate Saturated model to estimate means and variances
#      Two matrix styles - Two data styles
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


univSatModel1 <- mxModel("univSat1",
    manifestVars=selVars,
    mxPath(
        from=c("X"), 
        arrows=2, 
        free=T, 
        values=1, 
        lbound=.01, 
        labels="vX"
    ),
    mxData(
        observed=var(testData), 
        type="cov", 
        numObs=1000 
    ),
    type="RAM"
)
univSatFit1 <- mxRun(univSatModel1)
EC1 <- mxEval(S, univSatFit1)
LL1 <- mxEval(objective, univSatFit1)
SL1 <- univSatFit1@output$SaturatedLikelihood
Chi1 <- LL1-SL1
# example 1: Saturated Model with Cov Matrices and Path-Style Input
# -----------------------------------------------------------------------------

univSatModel1m <- mxModel("univSat1m",
    manifestVars=selVars,
    mxPath(
        from=c("X"), 
        arrows=2, 
        free=T, 
        values=1, 
        lbound=.01, 
        labels="vX"
    ),
    mxPath(
        from="one", 
        to="X", 
        arrows=1, 
        free=T, 
        values=0, 
        labels="mX"
    ),
    mxData(
        observed=var(testData), 
        type="cov", 
        numObs=1000, 
        means=colMeans(testData)
    ),
    type="RAM"
)
univSatFit1m <- mxRun(univSatModel1m)
EM1m <- mxEval(M, univSatFit1m)
EC1m <- mxEval(S, univSatFit1m)
LL1m <- mxEval(objective,univSatFit1m);
SL1m <- univSatFit1m@output$SaturatedLikelihood
Chi1m <- LL1m-SL1m
# example 1m: Saturated Model with Cov Matrices & Means and Path-Style Input
# -----------------------------------------------------------------------------

univSatModel2 <- mxModel("univSat2",
    manifestVars= selVars,
    mxPath(
        from=c("X"), 
        arrows=2, 
        free=T, 
        values=1, 
        lbound=.01, 
        labels="vX"
    ),
    mxPath(
        from="one", 
        to="X", 
        arrows=1, 
        free=T, 
        values=0, 
        labels="mX"
    ),    
    mxData(
        observed=testData, 
        type="raw", 
    ),
    type="RAM"
)
univSatFit2 <- mxRun(univSatModel2)
EM2 <- mxEval(M, univSatFit2)
EC2 <- mxEval(S, univSatFit2)
LL2 <- mxEval(objective,univSatFit2);
# example 2: Saturated Model with Raw Data and Path input
# -----------------------------------------------------------------------------

univSatModel2s <- mxModel(univSatModel1,
    mxData(
        observed=testData, 
        type="raw"
    ),
    mxPath(
        from="one", 
        to="X", 
        arrows=1, 
        free=T, 
        values=0, 
        labels="mX"
    ),    
    name="univSat2s",
    type="RAM"
    )
univSatFit2s <- mxRun(univSatModel2s)
EM2s <- mxEval(M, univSatFit2s)
EC2s <- mxEval(S, univSatFit2s)
LL2s <- mxEval(objective,univSatFit2s);
# example 2s: Saturated Model with Raw Data and Path input built upon 
# Cov/Means version
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
# example 3m: Saturated Model with Cov Matrices & Means and Matrix-Style Input
# -----------------------------------------------------------------------------

univSatModel4 <- mxModel("univSat4",
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
        observed=testData, 
        type="raw", 
    ),
    mxFIMLObjective(
        covariance="expCov", 
        means="expMean",
        dimnames=selVars
    )
)
univSatFit4 <- mxRun(univSatModel4)
EM4 <- mxEval(expMean, univSatFit4)
EC4 <- mxEval(expCov, univSatFit4)
LL4 <- mxEval(objective, univSatFit4);
# examples 4: Saturated Model with Raw Data and Matrix-Style Input
# -----------------------------------------------------------------------------

# example Mx..1: Saturated Model with Cov Matrices
# -------------------------------------
Mx.EC1 <-  1.062112
Mx.LL1 <- -1.474434e-17

# example Mx..1m: Saturated Model with Cov Matrices & Means
# -------------------------------------
Mx.EM1m <- 0.01680509
Mx.EC1m <- 1.062112
Mx.LL1m <- -1.108815e-13

Mx.EM2 <- 0.01680516
Mx.EC2 <- 1.061050
Mx.LL2 <- 2897.135
# Mx answers hard-coded
# -----------------------------------------------------------------------------

cov <- rbind(cbind(EC1,EC1m,EC2),cbind(EC3,EC3m,EC4))
mean <- rbind(cbind(EM1m, EM2),cbind(EM3m,EM4))
like <- rbind(cbind(LL1,LL1m,LL2),cbind(LL3,LL3m,LL4))
cov; mean; like
# OpenMx summary
# -----------------------------------------------------------------------------


Mx.cov <- cbind(Mx.EC1,Mx.EC1m,Mx.EC2)
Mx.mean <- cbind(Mx.EM1m,Mx.EM2)
Mx.like <- cbind(Mx.LL1,Mx.LL1m,Mx.LL2)
Mx.cov; Mx.mean; Mx.like
# old Mx summary
# -----------------------------------------------------------------------------

omxCheckCloseEnough(Chi1,Mx.LL1,.001)
omxCheckCloseEnough(EC1,Mx.EC1,.001)
# 1:CovPat
# -------------------------------------

omxCheckCloseEnough(Chi1m,Mx.LL1m,.001)
omxCheckCloseEnough(EC1m,Mx.EC1m,.001)
omxCheckCloseEnough(EM1m,Mx.EM1m,.001)
# 1m:CovMPat 
# -------------------------------------

omxCheckCloseEnough(LL2,Mx.LL2,.001)
omxCheckCloseEnough(EC2,Mx.EC2,.001)
omxCheckCloseEnough(EM2,Mx.EM2,.001)
# 2:RawPat 
# -------------------------------------

omxCheckCloseEnough(LL2s,Mx.LL2,.001)
omxCheckCloseEnough(EC2s,Mx.EC2,.001)
omxCheckCloseEnough(EM2s,Mx.EM2,.001)
# 2:RawSPat
# -------------------------------------

omxCheckCloseEnough(Chi3,Mx.LL1,.001)
omxCheckCloseEnough(EC3,Mx.EC1,.001)
# 3:CovMat
# -------------------------------------

omxCheckCloseEnough(Chi3m,Mx.LL1m,.001)
omxCheckCloseEnough(EC3m,Mx.EC1m,.001)
omxCheckCloseEnough(EM3m,Mx.EM1m,.001)
# 3m:CovMPat 
# -------------------------------------

omxCheckCloseEnough(LL4,Mx.LL2,.001)
omxCheckCloseEnough(EC4,Mx.EC2,.001)
omxCheckCloseEnough(EM4,Mx.EM2,.001)
# 4:RawMat
# -------------------------------------
# Compare OpenMx results to Mx results
# (LL: likelihood; EC: expected covariance, EM: expected means)
# -----------------------------------------------------------------------------
