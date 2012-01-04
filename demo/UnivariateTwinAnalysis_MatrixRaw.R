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
# Program: UnivariateTwinAnalysis_MatrixRaw.R  
# Author: Hermine Maes
# Date: 2009.08.01 
#
# ModelType: ACE
# DataType: Twin
# Field: Human Behavior Genetics
#
# Purpose: 
#      Univariate Twin Analysis model to estimate causes of variation
#      Matrix style model input - Raw data input
#
# RevisionHistory:
#      Hermine Maes -- 2009.10.08 updated & reformatted
#      Ross Gore -- 2011.06.06	added Model, Data & Field metadata
# -----------------------------------------------------------------------------

require(OpenMx)
# Load Library
# -----------------------------------------------------------------------------

data(twinData)
summary(twinData)
selVars <- c('bmi1','bmi2')
mzData <- as.matrix(subset(twinData, zyg==1, c(bmi1,bmi2)))
dzData <- as.matrix(subset(twinData, zyg==3, c(bmi1,bmi2)))
colMeans(mzData,na.rm=TRUE)
colMeans(dzData,na.rm=TRUE)
cov(mzData,use="complete")
cov(dzData,use="complete")
# Prepare Data
# -----------------------------------------------------------------------------

twinACEModel <- mxModel("twinACE",

	mxMatrix(
		type="Full", 
		nrow=1, 
		ncol=1, 
		free=TRUE,  
		values=.6,  
		label="a", 
		name="X"
	), 
	mxMatrix(
		type="Full", 
		nrow=1, 
		ncol=1, 
		free=TRUE,  
		values=.6,  
		label="c", 
		name="Y"
	),
	mxMatrix(
		type="Full", 
		nrow=1, 
		ncol=1, 
		free=TRUE,  
		values=.6,  
		label="e", 
		name="Z"
	),
	# Matrices X, Y, and Z to store 
	# a, c, and e path coefficients
	# ----------------------------------
			
	mxAlgebra(
		expression=X %*% t(X), 
		name="A"
	), 
	mxAlgebra(
		expression=Y %*% t(Y), 
		name="C"
	), 
	mxAlgebra(
		expression=Z %*% t(Z), 
		name="E"
	),
	mxMatrix(
		type="Full", 
		nrow=1, 
		ncol=2, 
		free=TRUE, 
		values= 20,
		label="mean", 
		name="expMean"
	),
	# Matrices A, C, and E compute 
	# variance components 
	# ----------------------------------
		
    mxAlgebra(
		expression= rbind  (cbind(A+C+E , A+C),
							cbind(A+C   , A+C+E)), 
		name="expCovMZ"
	),
	# Algebra for expected 
	# variance/covariance matrix in MZ
	# ----------------------------------
		
    mxAlgebra(
		expression= rbind  (cbind(A+C+E     , 0.5%x%A+C),
							cbind(0.5%x%A+C , A+C+E)), 
		name="expCovDZ"
	),
	mxModel("MZ",
	    mxData(
			observed=mzData, 
			type="raw"
		), 
	    mxFIMLObjective(
			covariance="twinACE.expCovMZ", 
			means="twinACE.expMean", 
			dimnames=selVars
		)
	),
	mxModel("DZ",
	    mxData(
			observed=dzData, 
			type="raw"
		), 
	    mxFIMLObjective(
			covariance="twinACE.expCovDZ", 
			means="twinACE.expMean", 
			dimnames=selVars
		)
	),
    mxAlgebra(
		expression=MZ.objective + DZ.objective, 
		name="twin"
	), 
	# Algebra for expected
	# variance/covariance matrix in DZ
	# note use of 0.5, converted 
	# to 1*1 matrix
	# ----------------------------------
	
	mxAlgebraObjective("twin")
)
# Fit ACE Model with RawData and Matrices Input
# -----------------------------------------------------------------------------

twinACEFit <- mxRun(twinACEModel)
# Run ACE model
# -----------------------------------------------------------------------------

MZc <- mxEval(expCovMZ, twinACEFit)  			# expected covariance matrix for MZ's
DZc <- mxEval(expCovDZ, twinACEFit)  			# expected covariance matrix for DZ's
M <- mxEval(expMean, twinACEFit)				# expected mean
A <- mxEval(a*a, twinACEFit)					# additive genetic variance, a^2
C <- mxEval(c*c, twinACEFit)					# shared environmental variance, c^2
E <- mxEval(e*e, twinACEFit)					# unique environmental variance, e^2
V <- (A+C+E)									# total variance
a2 <- A/V										# standardized additive genetic variance
c2 <- C/V										# standardized shared environmental variance
e2 <- E/V										# standardized unique environmental variance
ACEest <- rbind(cbind(A,C,E),cbind(a2,c2,e2)) 	# table of estimates
LL_ACE <- mxEval(objective, twinACEFit)			# likelihood of ACE model



Mx.A <- 0.6173023
Mx.C <- 5.595822e-14
Mx.E <- 0.1730462
Mx.M <- matrix(c(21.39293, 21.39293),1,2)
Mx.LL_ACE <- 4067.663
# 1: Heterogeneity Model
# -------------------------------------
# Mx answers hard-coded
# -----------------------------------------------------------------------------


omxCheckCloseEnough(LL_ACE,Mx.LL_ACE,.001)
omxCheckCloseEnough(A,Mx.A,.001)
omxCheckCloseEnough(C,Mx.C,.001)
omxCheckCloseEnough(E,Mx.E,.001)
omxCheckCloseEnough(M,Mx.M,.001)
# Compare OpenMx results to Mx results
# (LL: likelihood; EC: expected covariance, EM: expected means)
# -----------------------------------------------------------------------------


twinAEModel <- mxModel(twinACEModel, #twinAE model so we drop c at 0
	mxMatrix(
		type="Full", 
		nrow=1, 
		ncol=1, 
		free=FALSE, 
		values=0, 
		label="c", 
		name="Y"
	)	
)

twinAEFit <- mxRun(twinAEModel)
# Run AE model
# -----------------------------------------------------------------------------

MZc <- mxEval(expCovMZ, twinAEFit)
DZc <- mxEval(expCovDZ, twinAEFit)
M <- mxEval(expMean, twinAEFit)
A <- mxEval(a*a, twinAEFit)
C <- mxEval(c*c, twinAEFit)
E <- mxEval(e*e, twinAEFit)
V <- (A+C+E)
a2 <- A/V
c2 <- C/V
e2 <- E/V
AEest <- rbind(cbind(A,C,E),cbind(a2,c2,e2))
LL_AE <- mxEval(objective, twinAEFit)



Mx.A <- 0.6173023
Mx.C <- 0
Mx.E <- 0.1730462
Mx.M <- matrix(c(21.39293, 21.39293),1,2)
Mx.LL_AE <- 4067.663
# 1: Homogeneity Model
# -------------------------------------
# Hard-code Mx results
# -----------------------------------------------------------------------------

omxCheckCloseEnough(LL_AE,Mx.LL_AE,.001)
omxCheckCloseEnough(A,Mx.A,.001)
omxCheckCloseEnough(C,Mx.C,.001)
omxCheckCloseEnough(E,Mx.E,.001)
omxCheckCloseEnough(M,Mx.M,.001)

LRT_ACE_AE <- LL_AE - LL_ACE
# Compare OpenMx results to Mx results: 
# (LL: likelihood; EC: expected covariance, EM: expected means)
# -----------------------------------------------------------------------------

ACEest
AEest
LRT_ACE_AE
# Print relevant output
# -----------------------------------------------------------------------------