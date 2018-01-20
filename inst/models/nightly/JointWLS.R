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


#------------------------------------------------------------------------------
# Author: Michael D. Hunter
# Date: 2014-12-03
# Filename: JointWLS.R
# Purpose: Test WLS on Joint ordinal-continuous data.
#------------------------------------------------------------------------------



#------------------------------------------------
# a useful function

rms <- function(x, y=NA){
	if(is.matrix(x) && is.vector(y) && nrow(x) == length(y)){
		sqrt(colMeans((x-y)^2))
	} else if(is.matrix(x) && length(y) == 1 && is.na(y)){
		apply(x, 2, FUN=rms, x=x)
	} else{
		sqrt(mean((x-y)^2))
	}
}


#------------------------------------------------
# Load packages
# Generate data

require(mvtnorm)
require(Matrix)
require(OpenMx)

nvar <- 3
k    <- .5
sigma <- matrix(k, nvar, nvar)
diag(sigma) <- 1
N <- 20000 #500 # 20000

set.seed(943)


cDat    <- rmvnorm(N, sigma=sigma)
rawData <- (cDat>0) + (cDat>1)
rawData[rawData[,1]>1,1] <- 1

cCor <- cor(cDat)
oCor <- cor(rawData)

rawData <- data.frame(
	mxFactor(as.data.frame(rawData[,1]),      0:1),
	mxFactor(as.data.frame(rawData[,2:nvar]), 0:2)
	)
names(rawData) <- letters[(27-nvar):26]
nam <- names(rawData)


#------------------------------------------------
# Create joint data
# Create WLS data
# Fit Joint ML model
# Fit various least squares models

jdat <- data.frame(rawData, scale(cDat+rnorm(prod(dim(cDat)))))
names(jdat) <- c(nam, paste(nam, 'Con', sep=''))
r <- mxDataWLS(jdat)
u <- mxDataWLS(jdat, type="ULS")
d <- mxDataWLS(jdat, type="DLS")

jam <- names(jdat)

nvar <- ncol(jdat)
nord <- 3
nFactors <- 1
nThresholds <- 2


jod <- mxModel("jointThresholdModel",
	mxMatrix("Full", nvar, nFactors, values=0.2, free=TRUE, lbound=-.99, ubound=.99, name="L"),
	mxMatrix("Unit", nvar, 1, name="vectorofOnes"),
	mxMatrix("Zero", 1, nvar, name="M"),
	mxAlgebra(vectorofOnes - (diag2vec(L %*% t(L))) , name="E"),
	mxAlgebra(L %*% t(L) + vec2diag(E), name="impliedCovs"),
	mxMatrix("Full", nrow=nThresholds, ncol=nord, values=c(.2, NA, rep(c(.2, .4), 2)), free=c(TRUE, FALSE, rep(TRUE, 4)), name="thresholdMatrix", dimnames = list(c(), nam)),
	mxFitFunctionML(),
	mxExpectationNormal(covariance="impliedCovs", means="M", dimnames = jam, thresholds="thresholdMatrix", threshnames=nam),
	mxData(observed=jdat, type='raw'))

#ML
jodr <- mxRun(jod)

# WLS
#jow <- mxModel(jod, name="jointThesholdModelWls", r, mxFitFunctionWLS(),
#	mxExpectationNormal(covariance="impliedCovs", dimnames = jam, thresholds="thresholdMatrix", threshnames=nam))

jow <- mxModel(jod, name="jointThesholdModelWls", r, mxFitFunctionWLS())

jowr <- mxRun(jow)

# DWLS
jodw <- mxModel(jow, name="jointThesholdModelDwls", d)

jodwr <- mxRun(jodw)

#ULS
jou <- mxModel(jow, name="jointThesholdModelUls", u)

jour <- mxRun(jou)


round(sres <- cbind(ML=omxGetParameters(jodr), WLS=omxGetParameters(jowr), DWLS=omxGetParameters(jodwr), ULS=omxGetParameters(jour)), 4)

rms(sres)

# non full WLS estimates are very close to ML
omxCheckTrue(all(rms(sres)[-2,-2] < 0.02))

# full WLS estimates are off due to model misspecification
omxCheckTrue(all(rms(sres)[,2] < 0.15))


#------------------------------------------------------------------------------
# alternative joint model
# Check the standardization for Joint

jodm <- mxModel("jointThresholdModel",
	mxMatrix("Full", nvar, nFactors, values=0.2, free=TRUE, lbound=-.99, ubound=.99, name="L"),
	mxMatrix("Unit", nvar, 1, name="vectorofOnes"),
	mxMatrix("Full", 1, nvar, name="M", free=TRUE),
	mxAlgebra(vectorofOnes - (diag2vec(L %*% t(L))) , name="E"),
	mxAlgebra(L %*% t(L) + vec2diag(E), name="impliedCovs"),
	mxMatrix("Full", nThresholds, nord, values=c(-1, NA, -1, 0, -2, -1), name="thresholdMatrix", free=FALSE, dimnames = list(c(), nam)),
	mxFitFunctionML(),
	mxExpectationNormal(covariance="impliedCovs", means="M", dimnames = jam, thresholds="thresholdMatrix", threshnames=nam),
	mxData(observed=jdat, type='raw'))

jodmr <- mxRun(jodm)

jowm <- mxModel(jodm, name="jointThesholdModelWls", r, mxFitFunctionWLS())
jowmr <- mxRun(jowm)

# DWLS
jodwm <- mxModel(jowm, name="jointThesholdModelDwls", d)
jodwmr <- mxRun(jodwm)

#ULS
joum <- mxModel(jowm, name="jointThesholdModelUls", u)
joumr <- mxRun(joum)


round( stres <- cbind(ML=coef(jodmr), WLS=coef(jowmr), DWLS=coef(jodwmr), ULS=coef(joumr)), 4)

# Note that when the fixed thresholds are set to have a distance much GREATER or
#  much LESS than 1.0, I believe this contradicts the Total Variance of 1
#  "constraint" in mxAlgebra "E".  This makes the WLS/DWLS/ULS estimator diverge
#  from the ML.


rms(stres)

# non full WLS estimates are very close to ML
omxCheckTrue(all(rms(stres)[-2,-2] < 0.02))

# full WLS estimates are off due to model misspecification
omxCheckTrue(all(rms(stres)[,2] < 0.15))


#------------------------------------------------------------------------------
# End

