#
#   Copyright 2007-2014 The OpenMx Project
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
#      Hermine Maes -- 2014.11.02 piecewise specification
# -----------------------------------------------------------------------------

# Load Library
# -----------------------------------------------------------------------------

    require(OpenMx)

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

    
    # Set Starting Values
    svMe      <- 20      # start value for means
    svPa      <- .6      # start value for path coefficients (sqrt(variance/#ofpaths))

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
    covDZ     <- mxAlgebra( expression=rbind( cbind(V, 0.5%x%A + C), 
                                              cbind(0.5%x%A + C , V)), name="expCovDZ" )

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
    AceModel  <- mxModel( "ACE", pars, modelMZ, modelDZ, fitML )

    # Run ACE model
    AceFit    <- mxRun(AceModel, intervals=T)
    AceSumm   <- summary(AceFit)
    AceSumm

# Fit ACE Model with RawData and Matrices Input
# -----------------------------------------------------------------------------


    # Generate ACE Model Output
    estMean   <- mxEval(expMean, AceFit$MZ)       # expected mean
    estCovMZ  <- mxEval(expCovMZ, AceFit$MZ)      # expected covariance matrix for MZ's
    estCovDZ  <- mxEval(expCovDZ, AceFit$DZ)      # expected covariance matrix for DZ's
    estVA     <- mxEval(a*a, AceFit)              # additive genetic variance, a^2
    estVC     <- mxEval(c*c, AceFit)              # dominance variance, d^2
    estVE     <- mxEval(e*e, AceFit)              # unique environmental variance, e^2
    estVP     <- (estVA+estVC+estVE)              # total variance
    estPropVA <- estVA/estVP                      # standardized additive genetic variance
    estPropVC <- estVC/estVP                      # standardized dominance variance
    estPropVE <- estVE/estVP                      # standardized unique environmental variance
    estACE    <- rbind(cbind(estVA,estVC,estVE),  # table of estimates
                       cbind(estPropVA,estPropVC,estPropVE))
    LL_ACE    <- mxEval(objective, AceFit)        # likelihood of ADE model
    ACEest <- rbind(cbind(estVA,estVC,estVE),cbind(estPropVA,estPropVC,estPropVE))   # table of estimates

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
omxCheckCloseEnough(estVA,Mx.A,.001)
omxCheckCloseEnough(estVC,Mx.C,.001)
omxCheckCloseEnough(estVE,Mx.E,.001)
omxCheckCloseEnough(estMean,Mx.M,.001)
# Compare OpenMx results to Mx results
# (LL: likelihood; EC: expected covariance, EM: expected means)
# -----------------------------------------------------------------------------


    # Run AE model
    AeModel   <- mxModel( AceFit, name="AE" )
    AeModel   <- omxSetParameters( AeModel, labels="c11", free=FALSE, values=0 )
    AeFit     <- mxRun(AeModel)

# Run AE model
# -----------------------------------------------------------------------------

    # Generate AE Model Output
    estVA     <- mxEval(a*a, AeFit)               # additive genetic variance, a^2
    estVE     <- mxEval(e*e, AeFit)               # unique environmental variance, e^2
    estVP     <- (estVA+estVE)                    # total variance
    estPropVA <- estVA/estVP                      # standardized additive genetic variance
    estPropVE <- estVE/estVP                      # standardized unique environmental variance
    estAE     <- rbind(cbind(estVA,estVE),        # table of estimates
                       cbind(estPropVA,estPropVE))
    LL_AE     <- mxEval(objective, AeFit)         # likelihood of AE model
    AEest <- rbind(cbind(estVA,estVC,estVE),cbind(estPropVA,estPropVC,estPropVE))


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
omxCheckCloseEnough(estVA,Mx.A,.001)
omxCheckCloseEnough(estVC,Mx.C,.001)
omxCheckCloseEnough(estVE,Mx.E,.001)
omxCheckCloseEnough(estMean,Mx.M,.001)

	LRT_ACE_AE <- LL_AE - LL_ACE
# Compare OpenMx results to Mx results: 
# (LL: likelihood; EC: expected covariance, EM: expected means)
# -----------------------------------------------------------------------------

ACEest
AEest
LRT_ACE_AE
# Print relevant output
# -----------------------------------------------------------------------------
