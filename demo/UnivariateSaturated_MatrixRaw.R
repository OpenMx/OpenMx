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
# Program: UnivariateSaturated_MatrixRaw.R  
# Author: Hermine Maes
# Date: 2009.08.01 
#
# ModelType: Saturated
# DataType: Simulated
# Field: None
#
# Purpose: 
#      Univariate Saturated model to estimate means and variances
#      Matrix style model input - Raw data input
#
# RevisionHistory:
#      Hermine Maes -- 2009.10.08 updated & reformatted
#      Ross Gore -- 2011.06.06	added Model, Data & Field metadata
# -----------------------------------------------------------------------

require(OpenMx)
# Load Library
# -----------------------------------------------------------------------


set.seed(100)
x <- rnorm (1000, 0, 1)
testData <- as.matrix(x)
selVars <- c("X")
dimnames(testData) <- list(NULL, selVars)
summary(testData)
colMeans(testData)
var(testData)
# Simulate Data
# -----------------------------------------------------------------------

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
# -----------------------------------------------------------------------


Mx.EM2 <- 0.01680516
Mx.EC2 <- 1.061050
Mx.LL2 <- 2897.135
# Mx answers hard-coded
# example Mx..1: Saturated Model with Raw Data
# -----------------------------------------------------------------------


omxCheckCloseEnough(LL4,Mx.LL2,.001)
omxCheckCloseEnough(EC4,Mx.EC2,.001)
omxCheckCloseEnough(EM4,Mx.EM2,.001)
# 4:RawMat
# -------------------------------------
# Compare OpenMx results to Mx results 
# (LL: likelihood; EC: expected covariance, EM: expected means)
# -----------------------------------------------------------------------
