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

library(OpenMx)
#This test script doesn't need to be run with all 3 optimizers, and it runs fine with CSOLNP:
if(mxOption(NULL,"Default optimizer")!="CSOLNP"){stop("SKIP")}

data(twinData)

# Select Variables for Analysis
vars      <- 'bmi'                     # list of variables names
nv        <- 1                         # number of variables
ntv       <- nv*2                      # number of total variables
selVars   <- paste(vars,c(rep(1,nv),rep(2,nv)),sep="")

# Select Covariates for Analysis
twinData[,'age'] <- twinData[,'age']/100
twinData  <- twinData[-which(is.na(twinData$age)),]
covVars   <- 'age'

# Select Data for Analysis
mzData    <- subset(twinData, zyg==1, c(selVars, covVars))
dzData    <- subset(twinData, zyg==3, c(selVars, covVars))

# Set Starting Values
svBe      <- 0.01                      # start value for regressions
svMe      <- 20                        # start value for means
svPa      <- .2                        # start value for path coefficient
svPe      <- .5                        # start value for path coefficient for e

# ----------------------------------------------------------------------------------------------------------------------
# PREPARE MODEL

# Create Matrices for Covariates and linear Regression Coefficients
defL      <- mxMatrix( type="Full", nrow=1, ncol=1, free=FALSE, labels=c("data.age"), name="defL" )
pathBl    <- mxMatrix( type="Full", nrow=1, ncol=1, free=TRUE, values=svBe, label="b11", name="bl" )

# Create Algebra for expected Mean Matrices
meanG     <- mxMatrix( type="Full", nrow=1, ncol=ntv, free=TRUE, values=svMe, labels="b0", name="meanG" )
expMean   <- mxAlgebra( expression= meanG + cbind(defL%*%bl,defL%*%bl), name="expMean" )

# Create Matrices for Variance Components
covA      <- mxMatrix( type="Symm", nrow=nv, ncol=nv, free=TRUE, values=svPa, label="A11", name="A" ) 
covC      <- mxMatrix( type="Symm", nrow=nv, ncol=nv, free=TRUE, values=svPa, label="C11", name="C" )
covE      <- mxMatrix( type="Symm", nrow=nv, ncol=nv, free=TRUE, values=svPa, label="E11", name="E" )

# Create Algebra for expected Variance/Covariance Matrices in MZ & DZ twins
covP      <- mxAlgebra( expression= A+C+E, name="V" )
covMZ     <- mxAlgebra( expression= A+C, name="cMZ" )
covDZ     <- mxAlgebra( expression= 0.5%x%A+ C, name="cDZ" )
expCovMZ  <- mxAlgebra( expression= rbind( cbind(V, cMZ), cbind(t(cMZ), V)), name="expCovMZ" )
expCovDZ  <- mxAlgebra( expression= rbind( cbind(V, cDZ), cbind(t(cDZ), V)), name="expCovDZ" )

# Create Data Objects for Multiple Groups
dataMZ    <- mxData( observed=mzData, type="raw" )
dataDZ    <- mxData( observed=dzData, type="raw" )

# Create Expectation Objects for Multiple Groups
expMZ     <- mxExpectationNormal( covariance="expCovMZ", means="expMean", dimnames=selVars )
expDZ     <- mxExpectationNormal( covariance="expCovDZ", means="expMean", dimnames=selVars )
funML     <- mxFitFunctionML()

# Create Model Objects for Multiple Groups
pars      <- list( pathBl, meanG, covA, covC, covE, covP )
defs      <- list( defL )
modelMZ   <- mxModel( pars, defs, expMean, covMZ, expCovMZ, dataMZ, expMZ, funML, name="MZ" )
modelDZ   <- mxModel( pars, defs, expMean, covDZ, expCovDZ, dataDZ, expDZ, funML, name="DZ" )
multi     <- mxFitFunctionMultigroup( c("MZ","DZ") )

# Create Algebra for Variance Components
rowVC     <- rep('VC',nv)
colVC     <- rep(c('A','C','E','SA','SC','SE'),each=nv)
estVC     <- mxAlgebra( expression=cbind(A,C,E,A/V,C/V,E/V), name="VC", dimnames=list(rowVC,colVC) )

# Create Confidence Interval Objects
ciACE     <- mxCI( "VC[1,1:3]" )

# Build Model with Confidence Intervals
modelACE  <- mxModel( "oneACEvca", pars, modelMZ, modelDZ, multi, estVC, ciACE )

# Run ACE Model
fitACE    <- mxRun( modelACE, intervals=F )

#AE
modelAE   <- omxSetParameters( fitACE, labels="C11", free=FALSE, values=0, name="oneAEvca" )
fitAE     <- mxRun( modelAE, intervals=F )

#CE
modelCE   <- omxSetParameters( modelACE, labels="A11", free=FALSE, values=0, name="oneCEvca" )
fitCE     <- mxRun( modelCE, intervals=F )

#E
modelE    <- omxSetParameters( fitAE, labels="A11", free=FALSE, values=0, name="oneEvca" )
fitE <- mxRun(modelE)

#Test omxAkaikeWeights():
omxAkaikeWeights(models=list(fitACE,fitAE,fitCE,fitE))

#Test mxModelAverage():
mxModelAverage(reference="C11",models=list(fitACE,fitAE,fitCE,fitE))

#Test mxModelAverage(), with free parameters only (which is faster), 
#for all permutations of arguments `include`, `SE`, `refAsBlock`, & `type`:
mxModelAverage(
	reference=c("b11","b0","A11","C11","E11"),models=list(fitACE,fitAE,fitCE,fitE),include="onlyFree",SE=NULL,refAsBlock=FALSE,type="AIC")
mxModelAverage(
	reference=c("b11","b0","A11","C11","E11"),models=list(fitACE,fitAE,fitCE,fitE),include="all",SE=NULL,refAsBlock=FALSE,type="AIC")
mxModelAverage(
	reference=c("b11","b0","A11","C11","E11"),models=list(fitACE,fitAE,fitCE,fitE),include="onlyFree",SE=FALSE,refAsBlock=FALSE,type="AIC")
mxModelAverage(
	reference=c("b11","b0","A11","C11","E11"),models=list(fitACE,fitAE,fitCE,fitE),include="all",SE=FALSE,refAsBlock=FALSE,type="AIC")
mxModelAverage(
	reference=c("b11","b0","A11","C11","E11"),models=list(fitACE,fitAE,fitCE,fitE),include="onlyFree",SE=TRUE,refAsBlock=FALSE,type="AIC")
mxModelAverage(
	reference=c("b11","b0","A11","C11","E11"),models=list(fitACE,fitAE,fitCE,fitE),include="all",SE=TRUE,refAsBlock=FALSE,type="AIC")
omxCheckError(mxModelAverage(
		reference=c("b11","b0","A11","C11","E11"),models=list(fitACE,fitAE,fitCE,fitE),include="onlyFree",SE=NULL,refAsBlock=TRUE,type="AIC"),
		"when refAsBlock=TRUE and include='onlyFree', no references may be fixed in any model")
mxModelAverage(
	reference=c("b11","b0","A11","C11","E11"),models=list(fitACE,fitAE,fitCE,fitE),include="all",SE=NULL,refAsBlock=TRUE,type="AIC")
omxCheckError(mxModelAverage(
	reference=c("b11","b0","A11","C11","E11"),models=list(fitACE,fitAE,fitCE,fitE),include="onlyFree",SE=FALSE,refAsBlock=TRUE,type="AIC"),
	"when refAsBlock=TRUE and include='onlyFree', no references may be fixed in any model")
mxModelAverage(
	reference=c("b11","b0","A11","C11","E11"),models=list(fitACE,fitAE,fitCE,fitE),include="all",SE=FALSE,refAsBlock=TRUE,type="AIC")
omxCheckError(mxModelAverage(
	reference=c("b11","b0","A11","C11","E11"),models=list(fitACE,fitAE,fitCE,fitE),include="onlyFree",SE=TRUE,refAsBlock=TRUE,type="AIC"),
	"when refAsBlock=TRUE and include='onlyFree', no references may be fixed in any model")
mxModelAverage(
	reference=c("b11","b0","A11","C11","E11"),models=list(fitACE,fitAE,fitCE,fitE),include="all",SE=TRUE,refAsBlock=TRUE,type="AIC")
mxModelAverage(
	reference=c("b11","b0","A11","C11","E11"),models=list(fitACE,fitAE,fitCE,fitE),include="onlyFree",SE=NULL,refAsBlock=FALSE,type="AICc")
mxModelAverage(
	reference=c("b11","b0","A11","C11","E11"),models=list(fitACE,fitAE,fitCE,fitE),include="all",SE=NULL,refAsBlock=FALSE,type="AICc")
mxModelAverage(
	reference=c("b11","b0","A11","C11","E11"),models=list(fitACE,fitAE,fitCE,fitE),include="onlyFree",SE=FALSE,refAsBlock=FALSE,type="AICc")
mxModelAverage(
	reference=c("b11","b0","A11","C11","E11"),models=list(fitACE,fitAE,fitCE,fitE),include="all",SE=FALSE,refAsBlock=FALSE,type="AICc")
mxModelAverage(
	reference=c("b11","b0","A11","C11","E11"),models=list(fitACE,fitAE,fitCE,fitE),include="onlyFree",SE=TRUE,refAsBlock=FALSE,type="AICc")
mxModelAverage(
	reference=c("b11","b0","A11","C11","E11"),models=list(fitACE,fitAE,fitCE,fitE),include="all",SE=TRUE,refAsBlock=FALSE,type="AICc")
omxCheckError(mxModelAverage(
	reference=c("b11","b0","A11","C11","E11"),models=list(fitACE,fitAE,fitCE,fitE),include="onlyFree",SE=NULL,refAsBlock=TRUE,type="AICc"),
	"when refAsBlock=TRUE and include='onlyFree', no references may be fixed in any model")
mxModelAverage(
	reference=c("b11","b0","A11","C11","E11"),models=list(fitACE,fitAE,fitCE,fitE),include="all",SE=NULL,refAsBlock=TRUE,type="AICc")
omxCheckError(mxModelAverage(
	reference=c("b11","b0","A11","C11","E11"),models=list(fitACE,fitAE,fitCE,fitE),include="onlyFree",SE=FALSE,refAsBlock=TRUE,type="AICc"),
	"when refAsBlock=TRUE and include='onlyFree', no references may be fixed in any model")
mxModelAverage(
	reference=c("b11","b0","A11","C11","E11"),models=list(fitACE,fitAE,fitCE,fitE),include="all",SE=FALSE,refAsBlock=TRUE,type="AICc")
omxCheckError(mxModelAverage(
	reference=c("b11","b0","A11","C11","E11"),models=list(fitACE,fitAE,fitCE,fitE),include="onlyFree",SE=TRUE,refAsBlock=TRUE,type="AICc"),
	"when refAsBlock=TRUE and include='onlyFree', no references may be fixed in any model")
mxModelAverage(
	reference=c("b11","b0","A11","C11","E11"),models=list(fitACE,fitAE,fitCE,fitE),include="all",SE=TRUE,refAsBlock=TRUE,type="AICc")

#Test with MxAlgebra elements:
mxModelAverage(reference=c("V","MZ.expCovMZ[2,1]"),models=list(fitACE,fitAE,fitCE,fitE),include="onlyFree",SE=NULL,refAsBlock=FALSE,type="AIC")
mxModelAverage(reference=c("V","MZ.expCovMZ[2,1]"),models=list(fitACE,fitAE,fitCE,fitE),include="all",SE=NULL,refAsBlock=FALSE,type="AIC")
mxModelAverage(reference=c("V","MZ.expCovMZ[2,1]"),models=list(fitACE,fitAE,fitCE,fitE),include="onlyFree",SE=FALSE,refAsBlock=FALSE,type="AIC")
mxModelAverage(reference=c("V","MZ.expCovMZ[2,1]"),models=list(fitACE,fitAE,fitCE,fitE),include="all",SE=FALSE,refAsBlock=FALSE,type="AIC")
mxModelAverage(reference=c("V","MZ.expCovMZ[2,1]"),models=list(fitACE,fitAE,fitCE,fitE),include="onlyFree",SE=TRUE,refAsBlock=FALSE,type="AIC")
mxModelAverage(reference=c("V","MZ.expCovMZ[2,1]"),models=list(fitACE,fitAE,fitCE,fitE),include="all",SE=TRUE,refAsBlock=FALSE,type="AIC")
omxCheckError(mxModelAverage(
	reference=c("V","MZ.expCovMZ[2,1]"),models=list(fitACE,fitAE,fitCE,fitE),include="onlyFree",SE=NULL,refAsBlock=TRUE,type="AIC"),
	"when refAsBlock=TRUE and include='onlyFree', no references may be fixed in any model")
mxModelAverage(reference=c("V","MZ.expCovMZ[2,1]"),models=list(fitACE,fitAE,fitCE,fitE),include="all",SE=NULL,refAsBlock=TRUE,type="AIC")
omxCheckError(mxModelAverage(
	reference=c("V","MZ.expCovMZ[2,1]"),models=list(fitACE,fitAE,fitCE,fitE),include="onlyFree",SE=FALSE,refAsBlock=TRUE,type="AIC"),
	"when refAsBlock=TRUE and include='onlyFree', no references may be fixed in any model")
mxModelAverage(reference=c("V","MZ.expCovMZ[2,1]"),models=list(fitACE,fitAE,fitCE,fitE),include="all",SE=FALSE,refAsBlock=TRUE,type="AIC")
omxCheckError(mxModelAverage(
	reference=c("V","MZ.expCovMZ[2,1]"),models=list(fitACE,fitAE,fitCE,fitE),include="onlyFree",SE=TRUE,refAsBlock=TRUE,type="AIC"),
	"when refAsBlock=TRUE and include='onlyFree', no references may be fixed in any model")
mxModelAverage(reference=c("V","MZ.expCovMZ[2,1]"),models=list(fitACE,fitAE,fitCE,fitE),include="all",SE=TRUE,refAsBlock=TRUE,type="AIC")
mxModelAverage(reference=c("V","MZ.expCovMZ[2,1]"),models=list(fitACE,fitAE,fitCE,fitE),include="onlyFree",SE=NULL,refAsBlock=FALSE,type="AICc")
mxModelAverage(reference=c("V","MZ.expCovMZ[2,1]"),models=list(fitACE,fitAE,fitCE,fitE),include="all",SE=NULL,refAsBlock=FALSE,type="AICc")
mxModelAverage(
	reference=c("V","MZ.expCovMZ[2,1]"),models=list(fitACE,fitAE,fitCE,fitE),include="onlyFree",SE=FALSE,refAsBlock=FALSE,type="AICc")
mxModelAverage(
	reference=c("V","MZ.expCovMZ[2,1]"),models=list(fitACE,fitAE,fitCE,fitE),include="all",SE=FALSE,refAsBlock=FALSE,type="AICc")
mxModelAverage(
	reference=c("V","MZ.expCovMZ[2,1]"),models=list(fitACE,fitAE,fitCE,fitE),include="onlyFree",SE=TRUE,refAsBlock=FALSE,type="AICc")
mxModelAverage(
	reference=c("V","MZ.expCovMZ[2,1]"),models=list(fitACE,fitAE,fitCE,fitE),include="all",SE=TRUE,refAsBlock=FALSE,type="AICc")
omxCheckError(mxModelAverage(
	reference=c("V","MZ.expCovMZ[2,1]"),models=list(fitACE,fitAE,fitCE,fitE),include="onlyFree",SE=NULL,refAsBlock=TRUE,type="AICc"),
	"when refAsBlock=TRUE and include='onlyFree', no references may be fixed in any model")
mxModelAverage(
	reference=c("V","MZ.expCovMZ[2,1]"),models=list(fitACE,fitAE,fitCE,fitE),include="all",SE=NULL,refAsBlock=TRUE,type="AICc")
omxCheckError(mxModelAverage(
	reference=c("V","MZ.expCovMZ[2,1]"),models=list(fitACE,fitAE,fitCE,fitE),include="onlyFree",SE=FALSE,refAsBlock=TRUE,type="AICc"),
	"when refAsBlock=TRUE and include='onlyFree', no references may be fixed in any model")
mxModelAverage(
	reference=c("V","MZ.expCovMZ[2,1]"),models=list(fitACE,fitAE,fitCE,fitE),include="all",SE=FALSE,refAsBlock=TRUE,type="AICc")
omxCheckError(mxModelAverage(
	reference=c("V","MZ.expCovMZ[2,1]"),models=list(fitACE,fitAE,fitCE,fitE),include="onlyFree",SE=TRUE,refAsBlock=TRUE,type="AICc"),
	"when refAsBlock=TRUE and include='onlyFree', no references may be fixed in any model")
mxModelAverage(
	reference=c("V","MZ.expCovMZ[2,1]"),models=list(fitACE,fitAE,fitCE,fitE),include="all",SE=TRUE,refAsBlock=TRUE,type="AICc")
