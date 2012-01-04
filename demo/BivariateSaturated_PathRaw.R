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
# Program: BivariateSaturated_PathRaw.R  
# Author: Hermine Maes
# Date: 2009.08.01 
#
# ModelType: Saturated
# DataType: Continuous
# Field: None
#
# Purpose:
#      Bivariate Saturated model to estimate means and (co)variances
#      Path style model input - Raw data input
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
cov(testData)
# Simulate Data
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
    mxPath(
    	from="one",
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
LL2 <- mxEval(objective, bivSatFit2)
SL2 <- summary(bivSatFit2)$SaturatedLikelihood
Chi2 <- LL2-SL2
# Example 2: Saturated Model with Raw Data and Path input
# -----------------------------------------------------------------------------

Mx.EM2 <- matrix(c(0.03211188, -0.004889211),1,2)
Mx.EC2 <- matrix(c(1.0092891, 0.4813504, 0.4813504, 0.9935366),2,2)
Mx.LL2 <- 5415.772
# example Mx..2: Saturated Model with Raw Data
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
