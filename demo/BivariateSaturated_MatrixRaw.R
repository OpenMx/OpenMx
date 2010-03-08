#
#   Copyright 2007-2010 The OpenMx Project
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
# Program: BivariateSaturated_MatrixRaw.R  
#  Author: Hermine Maes
#    Date: 08 01 2009 
#
# Bivariate Saturated model to estimate means and (co)variances
# Matrix style model input - Raw data input
#
# Revision History
#   Hermine Maes -- 10 08 2009 updated & reformatted
# -----------------------------------------------------------------------

require(OpenMx)

#Simulate Data
# -----------------------------------------------------------------------
require(MASS)
set.seed(200)
rs=.5
xy <- mvrnorm (1000, c(0,0), matrix(c(1,rs,rs,1),2,2))
testData <- xy
selVars <- c("X","Y")
dimnames(testData) <- list(NULL, selVars)
summary(testData)
cov(testData)

#examples 4: Saturated Model with Raw Data and Matrix-Style Input
# -----------------------------------------------------------------------
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
	mxFIMLObjective(
	    covariance="expCov",
	    means="expMean",
	    dimnames=selVars
	)
)

bivSatFit4 <- mxRun(bivSatModel4)
bivSatSummary4 <- summary(bivSatFit4)
EM4 <- mxEval(expMean, bivSatFit4)
EC4 <- mxEval(expCov, bivSatFit4)
LL4 <- mxEval(objective,bivSatFit4)
omxCheckIdentical(bivSatSummary4$observedStatistics, sum(!is.na(testData)))

#Mx answers hard-coded
# -----------------------------------------------------------------------
#example Mx..2: Saturated Model with Raw Data
Mx.EM2 <- matrix(c(0.03211188, -0.004889211),1,2)
Mx.EC2 <- matrix(c(1.0092891, 0.4813504, 0.4813504, 0.9935366),2,2)
Mx.LL2 <- 5415.772


#Compare OpenMx results to Mx results 
# -----------------------------------------------------------------------
# (LL: likelihood; EC: expected covariance, EM: expected means)
#4:RawMat
omxCheckCloseEnough(LL4,Mx.LL2,.001)
omxCheckCloseEnough(EC4,Mx.EC2,.001)
omxCheckCloseEnough(EM4,Mx.EM2,.001)
