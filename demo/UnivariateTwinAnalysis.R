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
# Program: UnivariateTwinAnalysis.R  
# Author: Hermine Maes
# Date: 2009.08.01 
#
# ModelType: ACE
# DataType: Simulated Twin
# Field: Human Behavior Genetics
#
# Purpose: 
#      Univariate Twin Analysis Example in OpenMx:
#      From fitting saturated models to testing model assumptions 
#      To fitting the ACE model and a submodel
#
# RevisionHistory:
#      Hermine Maes -- 10 08 2009 updated & reformatted
#      Ross Gore -- 06 06 2011	added Model, Data & Field metadata
# -----------------------------------------------------------------------

require(OpenMx)
require(MASS)
# Load Library
# -----------------------------------------------------------------------

set.seed(200)
a2<-0.5		#Additive genetic variance component (a squared)
c2<-0.3		#Common environment variance component (c squared)
e2<-0.2		#Specific environment variance component (e squared)
rMZ <- a2+c2
rDZ <- .5*a2+c2
DataMZ <- mvrnorm (1000, c(0,0), matrix(c(1,rMZ,rMZ,1),2,2))
DataDZ <- mvrnorm (1000, c(0,0), matrix(c(1,rDZ,rDZ,1),2,2))

selVars <- c('t1','t2')
dimnames(DataMZ) <- list(NULL,selVars)
dimnames(DataDZ) <- list(NULL,selVars)
summary(DataMZ)
summary(DataDZ)
colMeans(DataMZ,na.rm=TRUE)
colMeans(DataDZ,na.rm=TRUE)
cov(DataMZ,use="complete")
cov(DataDZ,use="complete")
# Simulate Data: two standardized variables t1 & t2 for MZ's & DZ's
# -----------------------------------------------------------------------

twinSatModel <- mxModel("twinSat",
	mxModel("MZ",
		mxMatrix("Full", 1, 2, T, c(0,0), name="expMeanMZ"), 
		mxMatrix("Lower", 2, 2, T, .5, name="CholMZ"), 
		mxAlgebra(CholMZ %*% t(CholMZ), name="expCovMZ"), 
		mxData(DataMZ, type="raw"),
		mxFIMLObjective("expCovMZ", "expMeanMZ", selVars)),  
	mxModel("DZ",
		mxMatrix("Full", 1, 2, T, c(0,0), name="expMeanDZ"), 
		mxMatrix("Lower", 2, 2, T, .5, name="CholDZ"), 
		mxAlgebra(CholDZ %*% t(CholDZ), name="expCovDZ"), 
		mxData(DataDZ, type="raw"), 
		mxFIMLObjective("expCovDZ", "expMeanDZ", selVars)),
	mxAlgebra(MZ.objective + DZ.objective, name="twin"), 
	mxAlgebraObjective("twin"))
twinSatFit <- mxRun(twinSatModel, suppressWarnings=TRUE)
# Specify and Run Saturated Model with RawData and Matrix-style Input
# -----------------------------------------------------------------------


ExpMeanMZ <- mxEval(MZ.expMeanMZ, twinSatFit)
ExpCovMZ <- mxEval(MZ.expCovMZ, twinSatFit)
ExpMeanDZ <- mxEval(DZ.expMeanDZ, twinSatFit)
ExpCovDZ <- mxEval(DZ.expCovDZ, twinSatFit)
LL_Sat <- mxEval(objective, twinSatFit)
# Generate Saturated Model Output
# -----------------------------------------------------------------------

twinSatModelSub1 <- twinSatModel
twinSatModelSub1$MZ$expMeanMZ <- mxMatrix("Full", 1, 2, T, 0, "mMZ", name="expMeanMZ")
twinSatModelSub1$DZ$expMeanDZ <- mxMatrix("Full", 1, 2, T, 0, "mDZ", name="expMeanDZ")
twinSatFitSub1 <- mxRun(twinSatModelSub1, suppressWarnings=TRUE)
# Specify and Run Saturated SubModel 1 equating means across twin order
# -----------------------------------------------------------------------

twinSatModelSub2 <- twinSatModelSub1
twinSatModelSub2$MZ$expMeanMZ <- mxMatrix("Full", 1, 2, T, 0, "mean", name="expMeanMZ")
twinSatModelSub2$DZ$expMeanDZ <- mxMatrix("Full", 1, 2, T, 0, "mean", name="expMeanDZ")
twinSatFitSub2 <- mxRun(twinSatModelSub2, suppressWarnings=TRUE)
# Specify and Run Saturated SubModel 2 equating means across zygosity
# -----------------------------------------------------------------------

LL_Sat <- mxEval(objective, twinSatFit)
LL_Sub1 <- mxEval(objective, twinSatFitSub1)
LRT1 <- LL_Sub1 - LL_Sat
LL_Sub2 <- mxEval(objective, twinSatFitSub2)
LRT2 <- LL_Sub2 - LL_Sat
# Generate Saturated Model Comparison Output
# -----------------------------------------------------------------------

twinACEModel <- mxModel("twinACE", 
	mxMatrix("Full", 1, 2, T, 20, "mean", name="expMean"), 
	# Matrix expMean for expected mean 
	# vector for MZ and DZ twins    
	# -------------------------------------
	
	mxMatrix("Full", nrow=1, ncol=1, free=TRUE, values=.6, label="a", name="X"),
	mxMatrix("Full", nrow=1, ncol=1, free=TRUE, values=.6, label="c", name="Y"),
	mxMatrix("Full", nrow=1, ncol=1, free=TRUE, values=.6, label="e", name="Z"),
	# Matrices X, Y, and Z to store the 
	# a, c, and e path coefficients
	# -------------------------------------
	
	mxAlgebra(X * t(X), name="A"),
	mxAlgebra(Y * t(Y), name="C"),
	mxAlgebra(Z * t(Z), name="E"),	
	# Matrixes A, C, and E to compute 
	# A, C, and E variance components
	# -------------------------------------
	
	mxAlgebra(rbind(cbind(A+C+E   , A+C),
					cbind(A+C     , A+C+E)), name="expCovMZ"),
	# Matrix expCOVMZ for expected 
	# covariance matrix for MZ twins
	# -------------------------------------
	
	mxModel("MZ",
		mxData(DataMZ, type="raw"), 
		mxFIMLObjective("twinACE.expCovMZ", "twinACE.expMean",selVars)),

	mxAlgebra(rbind(cbind(A+C+E   , .5%x%A+C),
					cbind(.5%x%A+C , A+C+E)), name="expCovDZ"),
	# Matrix expCOVMZ for expected 
	# covariance matrix for DZ twins
	# -------------------------------------
	
	mxModel("DZ", 
		mxData(DataDZ, type="raw"), 
		mxFIMLObjective("twinACE.expCovDZ", "twinACE.expMean",selVars)),

	mxAlgebra(MZ.objective + DZ.objective, name="twin"), 
	mxAlgebraObjective("twin"))
twinACEFit <- mxRun(twinACEModel, suppressWarnings=TRUE)
# Specify and Run ACE Model with RawData and Matrix-style Input
# -----------------------------------------------------------------------

LL_ACE <- mxEval(objective, twinACEFit)
LRT_ACE= LL_ACE - LL_Sat
	MZc <- mxEval(expCovMZ, twinACEFit)
	DZc <- mxEval(expCovDZ, twinACEFit)
	M   <- mxEval(expMean, twinACEFit)
	# Retrieve expected mean vector and 
	# expected covariance matrices
	# -------------------------------------
	
	A <- mxEval(A, twinACEFit)
	C <- mxEval(C, twinACEFit)
	E <- mxEval(E, twinACEFit)
	# Retrieve the A, C, and E 
	# variance components
	# -------------------------------------
	
	V <- (A+C+E)
	a2 <- A/V
	c2 <- C/V
	e2 <- E/V
	# Calculate standardized variance 
	# components
	# -------------------------------------
	
	ACEest <- rbind(cbind(A,C,E),cbind(a2,c2,e2)) 
	ACEest <- data.frame(ACEest, row.names=c("Variance Components","Standardized VC"))
	names(ACEest)<-c("A", "C", "E")
 	ACEest; LL_ACE; LRT_ACE
	#Build and print reporting table with 
	# row and column names
	# -------------------------------------
	
# Generate ACE Model Output
# -----------------------------------------------------------------------
	
twinAEModel <- twinACEModel
twinAEModel$twinACE$Y <- mxMatrix("Full", 1, 1, F, 0, "c", name="Y")  # drop c
twinAEFit <- mxRun(twinAEModel, suppressWarnings=TRUE)
# Specify and Fit Reduced AE Model (Drop c @0)
# ----------------------------------------------------------------------

LL_AE <- mxEval(objective, twinAEFit)

	MZc <- mxEval(expCovMZ, twinAEFit)
	DZc <- mxEval(expCovDZ, twinAEFit)
	M   <- mxEval(expMean, twinAEFit)
	# Retrieve expected mean vector and 
	# expected covariance matrices
	# -------------------------------------
	A <- mxEval(A, twinAEFit)
	C <- mxEval(C, twinAEFit)
	E <- mxEval(E, twinAEFit)
	# Retrieve the A, C and E variance 
	# components
	# -------------------------------------
	
	V <- (A+C+E)
	a2 <- A/V
	c2 <- C/V
	e2 <- E/V
	# Calculate standardized variance 
	# components
	# -------------------------------------

	AEest <- rbind(cbind(A,C,E),cbind(a2,c2,e2)) 
	AEest <- data.frame(ACEest, row.names=c("Variance Components","Standardized VC"))
	names(ACEest)<-c("A", "C", "E")
	AEest; LL_AE; 
	# Build and print reporting table with 
	# row and column names
	# -------------------------------------

	LRT_ACE_AE <- LL_AE - LL_ACE
	LRT_ACE_AE
	# Calculate and print likelihood ratio 
	# test
	# -------------------------------------
# Generate ACE Model Output
# -----------------------------------------------------------------------