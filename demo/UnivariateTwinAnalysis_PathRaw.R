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
# Program: UnivariateTwinAnalysis_PathRaw.R  
#  Author: Hermine Maes
#    Date: 08 01 2009 
#
# Univariate Twin Analysis model to estimate causes of variation (ACE)
# Path style model input - Raw data input
#
# Revision History
#   Hermine Maes -- 10 08 2009 updated & reformatted
# -----------------------------------------------------------------------

require(OpenMx)

#Prepare Data
# -----------------------------------------------------------------------
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

#Fit ACE Model with RawData and Path-Style Input
# -----------------------------------------------------------------------
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

#Run ACE model
# -----------------------------------------------------------------------
twinACEFit <- mxRun(twinACEModel)

MZc <- mxEval(MZ.covariance, twinACEFit)
DZc <- mxEval(DZ.covariance, twinACEFit)
M <- mxEval(MZ.means, twinACEFit)
A <- mxEval(a*a, twinACEFit)
C <- mxEval(c*c, twinACEFit)
E <- mxEval(e*e, twinACEFit)
V <- (A+C+E)
a2 <- A/V
c2 <- C/V
e2 <- E/V
ACEest <- rbind(cbind(A,C,E),cbind(a2,c2,e2))
LL_ACE <- mxEval(objective, twinACEFit)


#Mx answers hard-coded
# -----------------------------------------------------------------------
#1: Heterogeneity Model
Mx.A <- 0.6173023
Mx.C <- 5.595822e-14
Mx.E <- 0.1730462
Mx.M <- matrix(c(21.39293, 21.39293),1,2)
Mx.LL_ACE <- 4067.663

#Compare OpenMx results to Mx results 
# -----------------------------------------------------------------------
# (LL: likelihood; EC: expected covariance, EM: expected means)
omxCheckCloseEnough(LL_ACE,Mx.LL_ACE,.001)
omxCheckCloseEnough(A,Mx.A,.001)
omxCheckCloseEnough(C,Mx.C,.001)
omxCheckCloseEnough(E,Mx.E,.001)
omxCheckCloseEnough(M,Mx.M,.001)


#Run AE model
# -----------------------------------------------------------------------
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

MZc <- mxEval(MZ.covariance, twinAEFit)
DZc <- mxEval(DZ.covariance, twinAEFit)
M <- mxEval(MZ.means, twinAEFit)
A <- mxEval(a*a, twinAEFit)
C <- mxEval(c*c, twinAEFit)
E <- mxEval(e*e, twinAEFit)
V <- (A+C+E)
a2 <- A/V
c2 <- C/V
e2 <- E/V
AEest <- rbind(cbind(A, C, E),cbind(a2, c2, e2))
LL_AE <- mxEval(objective, twinAEFit)

#Mx answers hard-coded
# -----------------------------------------------------------------------
#1: Homogeneity Model
Mx.A <- 0.6173023
Mx.C <- 0
Mx.E <- 0.1730462
Mx.M <- matrix(c(21.39293, 21.39293),1,2)
Mx.LL_AE <- 4067.663

#Compare OpenMx results to Mx results 
# -----------------------------------------------------------------------
# (LL: likelihood; EC: expected covariance, EM: expected means)
omxCheckCloseEnough(LL_AE,Mx.LL_AE,.001)
omxCheckCloseEnough(A,Mx.A,.001)
omxCheckCloseEnough(C,Mx.C,.001)
omxCheckCloseEnough(E,Mx.E,.001)
omxCheckCloseEnough(M,Mx.M,.001)

LRT_ACE_AE <- LL_AE - LL_ACE

#Print relevant output
# -----------------------------------------------------------------------
ACEest
AEest
LRT_ACE_AE
