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
# Program: UnivariateSaturated_PathCov.R  
# Author: Hermine Maes
# Date: 2009.08.01 
#
# ModelType: Saturated
# DataType: Simulated Continuous
# Field: None
#
# Purpose: 
#      Univariate Saturated model to estimate means and variances
#      Path style model input - Covariance matrix data input
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
    manifestVars= selVars,
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
SL1 <- summary(univSatFit1)$SaturatedLikelihood
Chi1 <- LL1 - SL1
# example 1: Saturated Model with Cov Matrices and Path-Style Input
# -----------------------------------------------------------------------------

univSatModel1m <- mxModel("univSat1m",
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
SL1m <- summary(univSatFit1m)$SaturatedLikelihood
Chi1m <- LL1m-SL1m
# example 1m: Saturated Model with Cov Matrices & Means and Path-Style input
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


# 1:CovPat
# -------------------------------------
omxCheckCloseEnough(Chi1,Mx.LL1,.001)
omxCheckCloseEnough(EC1,Mx.EC1,.001)
# 1m:CovMPat 
# -------------------------------------
omxCheckCloseEnough(Chi1m,Mx.LL1m,.001)
omxCheckCloseEnough(EC1m,Mx.EC1m,.001)
omxCheckCloseEnough(EM1m,Mx.EM1m,.001)
# Compare OpenMx results to Mx results
# (LL: likelihood; EC: expected covariance, EM: expected means) 
# -----------------------------------------------------------------------------
