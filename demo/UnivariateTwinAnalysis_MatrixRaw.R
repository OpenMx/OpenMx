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
#      Hermine Maes -- 2014.11.04 piecewise specification
# -----------------------------------------------------------------------------

require(OpenMx)
# Load Library
# -----------------------------------------------------------------------------

# Load Data
data(twinData)

# Select Variables for Analysis
Vars      <- 'bmi'
nv        <- 1       # number of variables
ntv       <- nv*2    # number of total variables
selVars   <- paste(Vars,c(rep(1,nv),rep(2,nv)),sep="")   #c('bmi1','bmi2')

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

require(OpenMx)

# Set Starting Values
svMe      <- 20      # start value for means
svPa      <- .5      # start value for path coefficients (sqrt(variance/#ofpaths))

# ACE Model
# Matrices declared to store a, d, and e Path Coefficients
pathA     <- mxMatrix( type="Full", nrow=nv, ncol=nv, 
                       free=TRUE, values=svPa, label="a11", name="a" ) 
pathC     <- mxMatrix( type="Full", nrow=nv, ncol=nv, 
                       free=TRUE, values=svPa, label="c11", name="c" )
pathE     <- mxMatrix( type="Full", nrow=nv, ncol=nv, 
                       free=TRUE, values=svPa, label="e11", name="e" )

# Matrices generated to hold A, C, and E computed Variance Components
covA      <- mxAlgebra( expression=a %*% t(a), name="A" )
covC      <- mxAlgebra( expression=c %*% t(c), name="C" ) 
covE      <- mxAlgebra( expression=e %*% t(e), name="E" )

# Algebra to compute total variances
covP      <- mxAlgebra( expression=A+C+E, name="V" )

# Algebra for expected Mean and Variance/Covariance Matrices in MZ & DZ twins
meanG     <- mxMatrix( type="Full", nrow=1, ncol=ntv, 
                       free=TRUE, values=svMe, label="mean", name="expMean" )
covMZ     <- mxAlgebra( expression=rbind( cbind(V, A+C), 
                                          cbind(A+C, V)), name="expCovMZ" )
covDZ     <- mxAlgebra( expression=rbind( cbind(V, 0.5%x%A+C), 
                                          cbind(0.5%x%A+C , V)), name="expCovDZ" )

# Data objects for Multiple Groups
dataMZ    <- mxData( observed=mzData, type="raw" )
dataDZ    <- mxData( observed=dzData, type="raw" )

# Objective objects for Multiple Groups
expMZ     <- mxExpectationNormal( covariance="expCovMZ", means="expMean", 
                                  dimnames=selVars )
expDZ     <- mxExpectationNormal( covariance="expCovDZ", means="expMean", 
                                  dimnames=selVars )
funML     <- mxFitFunctionML()

# Combine Groups
pars      <- list( pathA, pathC, pathE, covA, covC, covE, covP )
modelMZ   <- mxModel( pars, meanG, covMZ, dataMZ, expMZ, funML, name="MZ" )
modelDZ   <- mxModel( pars, meanG, covDZ, dataDZ, expDZ, funML, name="DZ" )
fitML     <- mxFitFunctionMultigroup(c("MZ.fitfunction","DZ.fitfunction") )
twinACEModel  <- mxModel( "ACE", pars, modelMZ, modelDZ, fitML )

# Run Model
twinACEFit   <- mxRun(twinACEModel, intervals=T)
twinACESum   <- summary(twinACEFit)
twinACESum
# Fit ACE Model with RawData and Matrices Input
# -----------------------------------------------------------------------------

twinACEFit <- mxRun(twinACEModel)
# Run ACE Model
# -----------------------------------------------------------------------------

# Generate ACE Model Output
estMean   <- mxEval(mean, twinACEFit)             # expected mean
estCovMZ  <- mxEval(expCovMZ, twinACEFit$MZ)      # expected covariance matrix for MZ's
estCovDZ  <- mxEval(expCovDZ, twinACEFit$DZ)      # expected covariance matrix for DZ's
estVA     <- mxEval(a*a, twinACEFit)              # additive genetic variance, a^2
estVC     <- mxEval(c*c, twinACEFit)              # dominance variance, d^2
estVE     <- mxEval(e*e, twinACEFit)              # unique environmental variance, e^2
estVP     <- (estVA+estVC+estVE)                  # total variance
estPropVA <- estVA/estVP                          # standardized additive genetic variance
estPropVC <- estVC/estVP                          # standardized dominance variance
estPropVE <- estVE/estVP                          # standardized unique environmental variance
estACE    <- rbind(cbind(estVA,estVC,estVE),      # table of estimates
                   cbind(estPropVA,estPropVC,estPropVE))
LL_ACE    <- mxEval(objective, twinACEFit)        # likelihood of ADE model
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
omxCheckCloseEnough(estVA,Mx.A,.001)
omxCheckCloseEnough(estVC,Mx.C,.001)
omxCheckCloseEnough(estVE,Mx.E,.001)
omxCheckCloseEnough(estMean,Mx.M,.001)
# Compare OpenMx results to Mx results
# -----------------------------------------------------------------------------


# Change Model
twinAEModel   <- mxModel( twinACEFit, name="AE" )
twinAEModel   <- omxSetParameters( twinAEModel, labels="c11", free=FALSE, values=0 )
twinAEFit     <- mxRun(twinAEModel)
# Run AE Model
# -----------------------------------------------------------------------------

# Generate AE Model Output
estVA     <- mxEval(a*a, twinAEFit)               # additive genetic variance, a^2
estVE     <- mxEval(e*e, twinAEFit)               # unique environmental variance, e^2
estVP     <- (estVA+estVE)                    # total variance
estPropVA <- estVA/estVP                      # standardized additive genetic variance
estPropVE <- estVE/estVP                      # standardized unique environmental variance
estAE     <- rbind(cbind(estVA,estVE),        # table of estimates
                   cbind(estPropVA,estPropVE))
LL_AE     <- mxEval(objective, twinAEFit)         # likelihood of AE model
LRT_ACE_AE <- LL_AE - LL_ACE
# Get Model Output
# -----------------------------------------------------------------------------

Mx.A <- 0.6173023
Mx.C <- 0
Mx.E <- 0.1730462
Mx.M <- 21.39293
Mx.LL_AE <- 4067.663
# Hard-code Mx results
# -----------------------------------------------------------------------------

omxCheckCloseEnough(LL_AE,Mx.LL_AE,.001)
omxCheckCloseEnough(estVA,Mx.A,.001)
omxCheckCloseEnough(estVC,Mx.C,.001)
omxCheckCloseEnough(estVE,Mx.E,.001)
omxCheckCloseEnough(estMean,Mx.M,.001)
# Compare OpenMx results to Mx results: 
# -----------------------------------------------------------------------------

estACE
estAE
LRT_ACE_AE
# Print relevant output
# -----------------------------------------------------------------------------
