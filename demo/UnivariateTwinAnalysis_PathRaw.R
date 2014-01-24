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
# Program: UnivariateTwinAnalysis_PathRaw.R  
# Author: Hermine Maes
# Date: 2009.08.01
#
# ModelType: ACE
# DataType: Twin
# Field: Human Behavior Genetics 
#
# Purpose: 
#      Univariate Twin Analysis model to estimate causes of variation 
#      Path style model input - Raw data input
#
# RevisionHistory:
#      Hermine Maes -- 2009.10.08 updated & reformatted
#      Ross Gore -- 2011.06.06	added Model, Data & Field metadata
# -----------------------------------------------------------------------------

require(OpenMx)
# Load Library
# -----------------------------------------------------------------------------

data(twinData)
twinVars <- c('fam','age','zyg','part','wt1','wt2','ht1','ht2','htwt1','htwt2','bmi1','bmi2')
#dimnames(twinData) <- list(NULL, twinVars)
summary(twinData)
selVars <- c('bmi1','bmi2')
aceVars <- c("A1","C1","E1","A2","C2","E2")
mzData <- as.matrix(subset(twinData, zyg==1, c(bmi1,bmi2)))
dzData <- as.matrix(subset(twinData, zyg==3, c(bmi1,bmi2)))
colMeans(mzData,na.rm=TRUE)
colMeans(dzData,na.rm=TRUE)
cov(mzData,use="complete")
cov(dzData,use="complete")
# Prepare Data
# -----------------------------------------------------------------------------

ACEModel <- mxModel("ACE", # Twin ACE Model -- Path Specification
	type="RAM",
    manifestVars=selVars,
    latentVars=aceVars,
    mxPath(
    	from=aceVars, 
    	arrows=2, 
    	free=FALSE, 
    	values=1
    ),
    mxPath(
    	from="one", 
    	to=aceVars, 
    	arrows=1, 
    	free=FALSE, 
    	values=0
    ),
    mxPath(
    	from="one", 
    	to=selVars, 
    	arrows=1, 
    	free=TRUE, 
    	values=20, 
    	labels= "mean"
    ),
    mxPath(
    	from=c("A1","C1","E1"), 
    	to="bmi1", 
    	arrows=1, 
    	free=TRUE, 
    	values=.6, 
    	label=c("a","c","e")
    ),
    mxPath(
    	from=c("A2","C2","E2"), 
    	to="bmi2", 
    	arrows=1, 
    	free=TRUE, 
    	values=.6, 
    	label=c("a","c","e")
    ),
    mxPath(
    	from="C1", 
    	to="C2", 
    	arrows=2, 
    	free=FALSE, 
    	values=1
    )
)    
mzModel <- mxModel(ACEModel, name="MZ",
	mxPath(
		from="A1", 
		to="A2", 
		arrows=2, 
		free=FALSE, 
		values=1
	),
	mxData(
		observed=mzData, 
		type="raw"
	)
)
dzModel <- mxModel(ACEModel, name="DZ", 
    mxPath(
    	from="A1", 
    	to="A2", 
    	arrows=2, 
    	free=FALSE, 
    	values=.5
    ),
    mxData(
    	observed=dzData, 
    	type="raw"
    )
)

twinACEModel <- mxModel("twinACE", mzModel, dzModel,
    mxAlgebra(
    	expression=MZ.objective + DZ.objective, 
    	name="twin"
    ), 
    mxAlgebraObjective("twin")
)
# Fit ACE Model with RawData and Path-Style Input
# -----------------------------------------------------------------------------

twinACEFit <- mxRun(twinACEModel)
# Run ACE model
# -----------------------------------------------------------------------------

MZc <- twinACEFit$MZ.objective@info$expCov
DZc <- twinACEFit$DZ.objective@info$expCov
M <- twinACEFit$MZ.objective@info$expMean
A <- mxEval(a*a, twinACEFit)
C <- mxEval(c*c, twinACEFit)
E <- mxEval(e*e, twinACEFit)
V <- (A+C+E)
a2 <- A/V
c2 <- C/V
e2 <- E/V
ACEest <- rbind(cbind(A,C,E),cbind(a2,c2,e2))
LL_ACE <- mxEval(objective, twinACEFit)



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
#Compare OpenMx results to Mx results 
# -----------------------------------------------------------------------------
# (LL: likelihood; EC: expected covariance, EM: expected means)


AEModel <- mxModel(ACEModel, #name="twinAE",
    mxPath(
    	from=c("A1","C1","E1"), 
    	to="bmi1", 
    	arrows=1, 
    	free=c(T,F,T),
    	values=c(.6,0,.6), 
    	label=c("a","c","e")
    ),
    mxPath(
    	from=c("A2","C2","E2"), 
    	to="bmi2", 
    	arrows=1, 
    	free=c(T,F,T),
    	values=c(.6,0,.6), 
    	label=c("a","c","e")
    )
)
mzModel <- mxModel(AEModel, name="MZ",
	mxPath(
		from="A1", 
		to="A2", 
		arrows=2, 
		free=FALSE, 
		values=1
	),
	mxData(
		observed=mzData, 
		type="raw"
	)
)
dzModel <- mxModel(AEModel, name="DZ", 
    mxPath(
    	from="A1", 
    	to="A2", 
    	arrows=2, 
    	free=FALSE, 
    	values=.5
    ),
    mxData(
    	observed=dzData, 
    	type="raw"
    )
)    
twinAEModel <- mxModel("twinAE", mzModel, dzModel,
    mxAlgebra(
    	expression=MZ.objective + DZ.objective, 
    	name="twin"
    ), 
    mxAlgebraObjective("twin")
)

twinAEFit <- mxRun(twinAEModel)
#Run AE model
# -----------------------------------------------------------------------------

MZc <- twinAEFit$MZ.objective@info$expCov
DZc <- twinAEFit$DZ.objective@info$expCov
M <- twinAEFit$MZ.objective@info$expMean
A <- mxEval(a*a, twinAEFit)
C <- mxEval(c*c, twinAEFit)
E <- mxEval(e*e, twinAEFit)
V <- (A+C+E)
a2 <- A/V
c2 <- C/V
e2 <- E/V
AEest <- rbind(cbind(A, C, E),cbind(a2, c2, e2))
LL_AE <- mxEval(objective, twinAEFit)



Mx.A <- 0.6173023
Mx.C <- 0
Mx.E <- 0.1730462
Mx.M <- matrix(c(21.39293, 21.39293),1,2)
Mx.LL_AE <- 4067.663
# 1: Homogeneity Model
# -------------------------------------
# Mx answers hard-coded
# -----------------------------------------------------------------------------

omxCheckCloseEnough(LL_AE,Mx.LL_AE,.001)
omxCheckCloseEnough(A,Mx.A,.001)
omxCheckCloseEnough(C,Mx.C,.001)
omxCheckCloseEnough(E,Mx.E,.001)
omxCheckCloseEnough(M,Mx.M,.001)
# Compare OpenMx results to Mx results 
# (LL: likelihood; EC: expected covariance, EM: expected means)
# -----------------------------------------------------------------------------


LRT_ACE_AE <- LL_AE - LL_ACE


ACEest
AEest
LRT_ACE_AE
#Print relevant output
# -----------------------------------------------------------------------------
