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
# Program: UnivariateSaturated_PathRaw.R  
# Author: Hermine Maes
# Date: 2009.08.01 
#
# ModelType: Saturated
# DataType: Simulated Continuous
# Field: None
#
# Purpose: 
#      Univariate Saturated model to estimate means and variances
#      Path style model input - Raw data input
#
# RevisionHistory:
#      Hermine Maes -- 2009.10.08 updated & reformatted
#      Ross Gore -- 2011.06.06	added Model, Data & Field
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
    	from=c("one"),
    	to=c("X"),
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
# example 2: Saturated Model with Raw Data and Path-Style Input
# -----------------------------------------------------------------------------



Mx.EM2 <- 0.01680516
Mx.EC2 <- 1.061050
Mx.LL2 <- 2897.135
# example Mx..2: Saturated Model 
# with Raw Data
# -------------------------------------
# Mx answers hard-coded
# -----------------------------------------------------------------------------



omxCheckCloseEnough(LL2,Mx.LL2,.001)
omxCheckCloseEnough(EC2,Mx.EC2,.001)
omxCheckCloseEnough(EM2,Mx.EM2,.001)
# 2:RawPat 
# -------------------------------------
# Compare OpenMx results to Mx results
# (LL: likelihood; EC: expected covariance, EM: expected means)
# -----------------------------------------------------------------------------
