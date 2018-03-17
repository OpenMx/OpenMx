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

require(OpenMx)
data(twinData)

Vars      <- 'bmi'
nv        <- 1       # number of variables
ntv       <- nv*2    # number of total variables
selVars   <- paste(Vars,c(rep(1,nv),rep(2,nv)),sep="")   #c('bmi1','bmi2')

# Select Data for Analysis
selData <- subset(twinData, twinData$zyg %in% c(1,3))
selData$mzdef <- as.numeric(selData$zyg==1)
selData$dzdef <- as.numeric(selData$zyg==3)
selData <- selData[,c(selVars,"mzdef","dzdef")]

mzData    <- subset(twinData, zyg==1, selVars)
dzData    <- subset(twinData, zyg==3, selVars)

# Generate Descriptive Statistics
colMeans(mzData,na.rm=TRUE)
colMeans(dzData,na.rm=TRUE)
cov(mzData,use="complete")
cov(dzData,use="complete")

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

# Use a different mean for MZ and DZ to test imxRowGradients() correctly handling different parameter labels in different submodels:
meanMZ     <- mxMatrix( type="Full", nrow=1, ncol=ntv, 
												free=TRUE, values=svMe, label="mzmean", name="MZmean" )
meanDZ     <- mxMatrix( type="Full", nrow=1, ncol=ntv, 
												free=TRUE, values=svMe, label="dzmean", name="DZmean" )
covMZ     <- mxAlgebra( expression=rbind( cbind(V, A+C), 
																					cbind(A+C, V)), name="expCovMZ" )
covDZ     <- mxAlgebra( expression=rbind( cbind(V, 0.5%x%A+C), 
																					cbind(0.5%x%A+C , V)), name="expCovDZ" )

# Data objects for Multiple Groups
dataMZ    <- mxData( observed=mzData, type="raw" )
dataDZ    <- mxData( observed=dzData, type="raw" )

# Objective objects for Multiple Groups
expMZ     <- mxExpectationNormal( covariance="expCovMZ", means="MZmean", 
																	dimnames=selVars )
expDZ     <- mxExpectationNormal( covariance="expCovDZ", means="DZmean", 
																	dimnames=selVars )
funML     <- mxFitFunctionML()

# Combine Groups
pars      <- list( pathA, pathC, pathE, covA, covC, covE, covP )
modelMZ   <- mxModel( pars, meanMZ, covMZ, dataMZ, expMZ, funML, name="MZ" )
modelDZ   <- mxModel( pars, meanDZ, covDZ, dataDZ, expDZ, funML, name="DZ" )
fitML     <- mxFitFunctionMultigroup(c("MZ.fitfunction","DZ.fitfunction") )
twinACEModel  <- mxModel( "ACE", pars, modelMZ, modelDZ, fitML )
twinACEFit <- mxRun(twinACEModel)
summary(twinACEFit)
multigroupRSE <- imxRobustSE(twinACEFit, details=TRUE)



#Single-group twin model:
singlegroup <- mxModel(
	"SingleGroupTwinModel",
	mxData(selData,type="raw"),
	pars, meanMZ, meanDZ,
	mxMatrix(type="Full",nrow=1,ncol=1,labels=c("data.mzdef"),name="MZdef"),
	mxMatrix(type="Full",nrow=1,ncol=1,labels=c("data.dzdef"),name="DZdef"),
	mxAlgebra(MZdef%x%MZmean + DZdef%x%DZmean, name="Mu"),
	mxAlgebra(MZdef + DZdef*0.5, name="kinship"),
	mxAlgebra(rbind(cbind(V, kinship%x%A+C), 
									cbind(kinship%x%A+C , V)), name="Sigma"),
	mxExpectationNormal(covariance="Sigma",means="Mu",dimnames=selVars),
	mxFitFunctionML()
)
singlegroupFit <- mxRun(singlegroup)
singlegroupRSE <- imxRobustSE(singlegroupFit,TRUE)
omxCheckCloseEnough(singlegroupRSE[[1]], multigroupRSE[[1]], 1e-5)
omxCheckCloseEnough(singlegroupRSE[[2]], multigroupRSE[[2]], 1e-6)
omxCheckCloseEnough(singlegroupRSE[[3]], multigroupRSE[[3]], 1e-5)
omxCheckCloseEnough(singlegroupRSE[[4]], multigroupRSE[[4]], 5e-3)
omxCheckCloseEnough(singlegroupRSE[[5]], multigroupRSE[[5]], 1e-4)



#Check for existence and dimnames of output:
pn <- c("a11","c11","e11","mzmean","dzmean")
omxCheckTrue(length(multigroupRSE$SE))
omxCheckEquals(names(multigroupRSE$SE),pn)
omxCheckTrue(length(multigroupRSE$cov))
omxCheckEquals(rownames(multigroupRSE$cov),pn)
omxCheckEquals(colnames(multigroupRSE$cov),pn)
omxCheckTrue(length(multigroupRSE$bread))
omxCheckEquals(rownames(multigroupRSE$bread),pn)
omxCheckEquals(colnames(multigroupRSE$bread),pn)
omxCheckTrue(length(multigroupRSE$meat))
omxCheckEquals(rownames(multigroupRSE$meat),pn)
omxCheckEquals(colnames(multigroupRSE$meat),pn)
omxCheckTrue(length(multigroupRSE$TIC))
omxCheckTrue(length(singlegroupRSE$SE))
omxCheckEquals(names(singlegroupRSE$SE),pn)
omxCheckTrue(length(singlegroupRSE$cov))
omxCheckEquals(rownames(singlegroupRSE$cov),pn)
omxCheckEquals(colnames(singlegroupRSE$cov),pn)
omxCheckTrue(length(singlegroupRSE$bread))
omxCheckEquals(rownames(singlegroupRSE$bread),pn)
omxCheckEquals(colnames(singlegroupRSE$bread),pn)
omxCheckTrue(length(singlegroupRSE$meat))
omxCheckEquals(rownames(singlegroupRSE$meat),pn)
omxCheckEquals(colnames(singlegroupRSE$meat),pn)
omxCheckTrue(length(singlegroupRSE$TIC))
