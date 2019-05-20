#
#   Copyright 2007-2018 by the individuals mentioned in the source code history
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

# This test is adapted from Hermine's oneACEo.R .
# It's a monophenotype ACE model with an ordinal trait, and uses an equality MxConstraint to identify the model. 
# The script deliberately puts the MxConstraint into the container MxModel, the MZ submodel, and the DZ submodel, 
# thereby making it redundant. 

library(OpenMx)
#library(polycor)

# ----------------------------------------------------------------------------------------------------------------------
# PREPARE DATA

# Load Data
data(twinData)
dim(twinData)

# Create Ordinal Variables
nth       <- 3                         # number of thresholds
ordData   <- data.frame(obmi1= twinData[,'bmi1'],
												obmi2= twinData[,'bmi2'], zyg= twinData[,'zyg'])
quant     <- quantile(ordData[,c('obmi1','obmi2')],(0:(nth+1))/(nth+1),na.rm=TRUE)
for (i in c('obmi1','obmi2')) { ordData[,i] <- cut(ordData[,i], breaks=quant, labels=c(0:nth)) }

# Select Variables for Analysis
vars      <- 'obmi'                    # list of variables names
nv        <- 1                         # number of variables
ntv       <- nv*2                      # number of total variables
selVars   <- paste(vars,c(rep(1,nv),rep(2,nv)),sep="")

# Select Data for Analysis
mzData    <- subset(ordData, zyg==1, selVars)
dzData    <- subset(ordData, zyg==3, selVars)
mzDataF   <- mxFactor( x=mzData, levels=c(0:nth) )
dzDataF   <- mxFactor( x=dzData, levels=c(0:nth) )

# Generate Descriptive Statistics
sapply(mzData,table)
sapply(dzData,table)

# Set Starting Values
svLTh     <- -1.5                      # start value for first threshold
svITh     <- 1                         # start value for increments
svTh      <- matrix(rep(c(svLTh,(rep(svITh,nth-1)))),nrow=nth,ncol=nv)         # start value for thresholds
lbTh      <- matrix(rep(c(-3,(rep(0.001,nth-1))),nv),nrow=nth,ncol=nv)         # lower bounds for thresholds
svPa      <- .5                        # start value for path coefficient
svPe      <- .7                        # start value for path coefficient for e
lbPa      <- .0001                     # lower bound for path coefficient

# ----------------------------------------------------------------------------------------------------------------------
# PREPARE MODEL

# Create Algebra for expected Mean & Threshold Matrices
meanG     <- mxMatrix( type="Zero", nrow=1, ncol=ntv, name="meanG" )
thinG     <- mxMatrix( 
	type="Full", nrow=nth, ncol=ntv, free=TRUE, values=svTh, lbound=lbTh, 
	labels=c("t11","t12","t13","t21","t22","t23"), name="thinG" )
inc       <- mxMatrix( type="Lower", nrow=nth, ncol=nth, free=FALSE, values=1, name="inc" )
threG     <- mxAlgebra( expression= inc %*% thinG, name="threG" )

# Create Matrices for Path Coefficients
pathA     <- mxMatrix( type="Lower", nrow=nv, ncol=nv, free=TRUE, values=svPa, label="a11", lbound=lbPa, name="a" ) 
pathC     <- mxMatrix( type="Lower", nrow=nv, ncol=nv, free=TRUE, values=svPa, label="c11", lbound=lbPa, name="c" )
pathE     <- mxMatrix( type="Lower", nrow=nv, ncol=nv, free=TRUE, values=svPe, label="e11", lbound=lbPa, name="e" )

# Create Algebra for Variance Components
covA      <- mxAlgebra( expression=a %*% t(a), name="A" )
covC      <- mxAlgebra( expression=c %*% t(c), name="C" ) 
covE      <- mxAlgebra( expression=e %*% t(e), name="E" )

# Create Algebra for expected Variance/Covariance Matrices in MZ & DZ twins
covP      <- mxAlgebra( expression= A+C+E, name="V" )
covMZ     <- mxAlgebra( expression= A+C, name="cMZ" )
covDZ     <- mxAlgebra( expression= 0.5%x%A+ C, name="cDZ" )
expCovMZ  <- mxAlgebra( expression= rbind( cbind(V, cMZ), cbind(t(cMZ), V)), name="expCovMZ" )
expCovDZ  <- mxAlgebra( expression= rbind( cbind(V, cDZ), cbind(t(cDZ), V)), name="expCovDZ" )

# Constrain Variance of Binary Variables
var1     <- mxConstraint( expression=diag2vec(V)==1, name="Var1" )

# Create Data Objects for Multiple Groups
dataMZ    <- mxData( observed=mzDataF, type="raw" )
dataDZ    <- mxData( observed=dzDataF, type="raw" )

# Create Expectation Objects for Multiple Groups
expMZ     <- mxExpectationNormal( covariance="expCovMZ", means="meanG", dimnames=selVars, thresholds="threG" )
expDZ     <- mxExpectationNormal( covariance="expCovDZ", means="meanG", dimnames=selVars, thresholds="threG" )
funML     <- mxFitFunctionML()

# Create Model Objects for Multiple Groups
pars      <- list( meanG, thinG, inc, threG, covA, pathA, pathC, pathE, covC, covE, covP )
modelMZ   <- mxModel( pars, covMZ, expCovMZ, dataMZ, expMZ, funML, var1, name="MZ" )
modelDZ   <- mxModel( pars, covDZ, expCovDZ, dataDZ, expDZ, funML, var1, name="DZ" )
multi     <- mxFitFunctionMultigroup( c("MZ","DZ") )

# Create Algebra for Variance Components
rowVC     <- rep('VC',nv)
colVC     <- rep(c('A','C','E','SA','SC','SE'),each=nv)
estVC     <- mxAlgebra( expression=cbind(A,C,E,A/V,C/V,E/V), name="VC", dimnames=list(rowVC,colVC) )

# Create Confidence Interval Objects
ciACE     <- mxCI( "VC[1,1:3]" )

# Build Model with Confidence Intervals
modelACE  <- mxModel( "oneACEo", pars, var1, modelMZ, modelDZ, multi, estVC, ciACE )

if(mxOption(NULL,"Default optimizer")!="CSOLNP"){
	fitACE    <- mxRun( modelACE, intervals=F )
	summary(fitACE)
	
	modelACE2 <- modelACE
	modelACE2$MZ@constraints <- list()
	modelACE2$DZ@constraints <- list()
	fitACE2 <- mxRun(modelACE2,intervals=F)
	summary(fitACE2)
	
	omxCheckCloseEnough(coef(fitACE), coef(fitACE2), 1e-7)
	omxCheckCloseEnough(fitACE$output$standardErrors, fitACE2$output$standardErrors, 1e-7)
	omxCheckCloseEnough(fitACE$output$fit, fitACE2$output$fit, 1e-7)
} else{
	plan <- omxDefaultComputePlan()
	plan$steps <- list(GD=plan$steps$GD)
	plan$steps$GD <- mxComputeNelderMead(eqConstraintMthd="GDsearch")
	modelACE <- mxModel(modelACE, plan)
	omxCheckWarning(
		mxRun(modelACE),
		"counted 3 equality constraints, but equality-constraint Jacobian is apparently rank 1 at the start values; Nelder-Mead will not work correctly unless equality constraints are linearly independent (this warning may be spurious if there are non-smooth equality constraints)"
	)
	#fitACE <- mxRun(modelACE)
	#summary(fitACE)
}
