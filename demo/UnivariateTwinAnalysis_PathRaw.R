#
#   Copyright 2007-2018 The OpenMx Project
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
#      Hermine Maes -- 2014.11.04 piecewise specification
# -----------------------------------------------------------------------------

require(OpenMx)
# Load Library
# -----------------------------------------------------------------------------

# Load Data
data(twinData)

# Select Variables for Analysis
selVars   <- c('bmi1','bmi2')
aceVars   <- c("A1","C1","E1","A2","C2","E2")

# Select Data for Analysis
mzData    <- subset(twinData, zyg==1, selVars)
dzData    <- subset(twinData, zyg==3, selVars)

# Generate Descriptive Statistics
colMeans(mzData,na.rm=TRUE)
colMeans(dzData,na.rm=TRUE)
cov(mzData,use="complete")
cov(dzData,use="complete")
# Prepare Data
# -----------------------------------------------------------------------------

# Path objects for Multiple Groups
manifestVars=selVars
latentVars=aceVars
# variances of latent variables
latVariances <- mxPath( from=aceVars, arrows=2, 
                        free=FALSE, values=1 )
# means of latent variables
latMeans     <- mxPath( from="one", to=aceVars, arrows=1, 
                        free=FALSE, values=0 )
# means of observed variables
obsMeans     <- mxPath( from="one", to=selVars, arrows=1, 
                        free=TRUE, values=20, labels="mean" )
# path coefficients for twin 1
pathAceT1    <- mxPath( from=c("A1","C1","E1"), to="bmi1", arrows=1, 
                        free=TRUE, values=.5,  label=c("a","c","e") )
# path coefficients for twin 2
pathAceT2    <- mxPath( from=c("A2","C2","E2"), to="bmi2", arrows=1, 
                        free=TRUE, values=.5,  label=c("a","c","e") )
# covariance between C1 & C2
covC1C2      <- mxPath( from="C1", to="C2", arrows=2, 
                        free=FALSE, values=1 )
# covariance between A1 & A2 in MZ twins
covA1A2_MZ   <- mxPath( from="A1", to="A2", arrows=2, 
                        free=FALSE, values=1 )
# covariance between A1 & A2 in DZ twins
covA1A2_DZ   <- mxPath( from="A1", to="A2", arrows=2, 
                        free=FALSE, values=.5 )

# Data objects for Multiple Groups
dataMZ       <- mxData( observed=mzData, type="raw" )
dataDZ       <- mxData( observed=dzData, type="raw" )

# Combine Groups
paths        <- list( latVariances, latMeans, obsMeans, 
                      pathAceT1, pathAceT2, covC1C2 )
modelMZ      <- mxModel(model="MZ", type="RAM", manifestVars=selVars, 
                        latentVars=aceVars, paths, covA1A2_MZ, dataMZ )
modelDZ      <- mxModel(model="DZ", type="RAM", manifestVars=selVars, 
                        latentVars=aceVars, paths, covA1A2_DZ, dataDZ )
minus2ll     <- mxAlgebra( expression=MZ.fitfunction + DZ.fitfunction, 
                           name="minus2loglikelihood" )
obj          <- mxFitFunctionAlgebra( "minus2loglikelihood" )
twinACEModel     <- mxModel(model="ACE", modelMZ, modelDZ, minus2ll, obj )

# Run Model
twinACEFit   <- mxRun(twinACEModel)
twinACESum   <- summary(twinACEFit)
twinACESum
# Fit ACE Model with RawData and Path-Style Input
# -----------------------------------------------------------------------------

# Generate & Print Output
# mean
M  <- mxEval(mean, twinACEFit)
# additive genetic variance, a^2
A  <- mxEval(a*a, twinACEFit)
# shared environmental variance, c^2
C  <- mxEval(c*c, twinACEFit)
# unique environmental variance, e^2
E  <- mxEval(e*e, twinACEFit)
# total variance
V  <- (A+C+E)
# standardized A
a2 <- A/V
# standardized C
c2 <- C/V
# standardized E
e2 <- E/V
# table of estimates
estACE <- rbind(cbind(A,C,E),cbind(a2,c2,e2))
# likelihood of ACE model
LL_ACE <- mxEval(fitfunction, twinACEFit)
# Get Model Output
# -----------------------------------------------------------------------------


Mx.A <- 0.6173023
Mx.C <- 5.595822e-14
Mx.E <- 0.1730462
Mx.M <- 21.39293
Mx.LL_ACE <- 4067.663
# Mx answers hard-coded
# -----------------------------------------------------------------------------


omxCheckCloseEnough(LL_ACE,Mx.LL_ACE,.001)
omxCheckCloseEnough(A,Mx.A,.001)
omxCheckCloseEnough(C,Mx.C,.001)
omxCheckCloseEnough(E,Mx.E,.001)
omxCheckCloseEnough(M,Mx.M,.001)
#Compare OpenMx results to Mx results 
# -----------------------------------------------------------------------------


# Change Path Objects
# path coefficients for twin 1
pathAceT1    <- mxPath( from=c("A1","C1","E1"), to="bmi1", arrows=1, 
                        free=c(T,F,T), values=c(.6,0,.6),  label=c("a","c","e") )
# path coefficients for twin 2
pathAceT2    <- mxPath( from=c("A2","C2","E2"), to="bmi2", arrows=1, 
                        free=c(T,F,T), values=c(.6,0,.6),  label=c("a","c","e") )

# Combine Groups
paths        <- list( latVariances, latMeans, obsMeans, 
                        pathAceT1, pathAceT2, covC1C2 )
modelMZ      <- mxModel(model="MZ", type="RAM", manifestVars=selVars, 
                        latentVars=aceVars, paths, covA1A2_MZ, dataMZ )
modelDZ      <- mxModel(model="DZ", type="RAM", manifestVars=selVars, 
                        latentVars=aceVars, paths, covA1A2_DZ, dataDZ )
twinAEModel      <- mxModel(model="AE", modelMZ, modelDZ, minus2ll, obj )

# Run Model
twinAEFit    <- mxRun(twinAEModel)
twinAESum    <- summary(twinAEFit)
# Fit AE model
# -----------------------------------------------------------------------------

# Generate & Print Output
M  <- mxEval(mean, twinAEFit)
A  <- mxEval(a*a, twinAEFit)
C  <- mxEval(c*c, twinAEFit)
E  <- mxEval(e*e, twinAEFit)
V  <- (A+C+E)
a2 <- A/V
c2 <- C/V
e2 <- E/V
estAE <- rbind(cbind(A, C, E),cbind(a2, c2, e2))
LL_AE <- mxEval(fitfunction, twinAEFit)
LRT_ACE_AE <- LL_AE - LL_ACE
estACE
estAE
LRT_ACE_AE
# Get Model Output
# -----------------------------------------------------------------------------

Mx.A <- 0.6173023
Mx.C <- 0
Mx.E <- 0.1730462
Mx.M <- 21.39293
Mx.LL_AE <- 4067.663
# Mx answers hard-coded
# -----------------------------------------------------------------------------

omxCheckCloseEnough(LL_AE,Mx.LL_AE,.001)
omxCheckCloseEnough(A,Mx.A,.001)
omxCheckCloseEnough(C,Mx.C,.001)
omxCheckCloseEnough(E,Mx.E,.001)
omxCheckCloseEnough(M,Mx.M,.001)
# Compare OpenMx results to Mx results 
# -----------------------------------------------------------------------------

estACE
estAE
LRT_ACE_AE
#Print relevant output
# -----------------------------------------------------------------------------

