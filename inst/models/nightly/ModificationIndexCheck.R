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
# Date: 2015-03-11
# Filename: ModificationIndexCheck.R
# Purpose: Check that modification indices are behaving reasonably.
#------------------------------------------------------------------------------



#------------------------------------------------------------------------------
# Factor model check

require(OpenMx)
data(demoOneFactor)
manifests <- names(demoOneFactor)
latents <- c("G")
factorModel <- mxModel("One Factor",
      type="RAM",
      manifestVars = manifests,
      latentVars = latents,
      mxPath(from=latents, to=manifests),
      mxPath(from=manifests, arrows=2),
      mxPath(from=latents, arrows=2,
            free=FALSE, values=1.0),
      mxPath(from = 'one', to = manifests),
      mxData(demoOneFactor, type="raw"))
factorRun <- mxRun(factorModel)

fim <- mxMI(factorRun)


mplusON.mi <- c(NA, 0.001, 1.582, .272, .285, NA, 0.001, NA, 1.540, 4.041, .748, NA, 1.582, 1.540, NA, .151, 5.339, NA, .272, 4.042, .151, NA, 3.574, NA, NA, .748, 5.339, 3.574, NA, NA, NA)
mplusWITH.mi <- c(0.001, 1.582, 1.540, .272, NA, .151, .285, .748, 5.339, 3.574)

mplus.mi <- c(mplusON.mi, mplusWITH.mi)

plot(fim$MI.Full[1:length(mplus.mi)], mplus.mi)
abline(a=0, b=1)
points(fim$MI[1:length(mplus.mi)], mplus.mi, col='blue')
(theCor <- cor(fim$MI.Full[1:length(mplus.mi)], mplus.mi, use="pair"))
rms <- function(x, y){sqrt(mean((x-y)^2, na.rm=TRUE))}
(theRMS <- rms(fim$MI.Full[1:length(mplus.mi)], mplus.mi))

omxCheckTrue(theCor > 0.95)
omxCheckTrue(theRMS < 0.50)

omCheck <- fim$MI.Full[1:length(mplus.mi)][!is.na(mplus.mi)]
mpCheck <- mplus.mi[!is.na(mplus.mi)]
relMean <- mean(omCheck/mpCheck)
omxCheckTrue( relMean < 1.05 & relMean > 0.95)


#------------------------------------------------------------------------------
# State space model check


# Create data based on state space model.
require(OpenMx)
nvar <- 5
varnames <- paste("x", 1:nvar, sep="")
ssModel <- mxModel(model="State Space Manual Example",
    mxMatrix("Full", 1, 1, TRUE, .3, name="A"),
    mxMatrix("Zero", 1, 1, name="B"),
    mxMatrix("Full", nvar, 1, TRUE, .6, name="C", dimnames=list(varnames, "F1")),
    mxMatrix("Zero", nvar, 1, name="D"),
    mxMatrix("Diag", 1, 1, FALSE, 1, name="Q"),
    mxMatrix("Diag", nvar, nvar, TRUE, .2, name="R"),
    mxMatrix("Zero", 1, 1, name="x0"),
    mxMatrix("Diag", 1, 1, FALSE, 1, name="P0"),
    mxMatrix("Zero", 1, 1, name="u"),
    mxExpectationStateSpace("A", "B", "C", "D", "Q", "R", "x0", "P0", "u"),
    mxFitFunctionML()
)


Rmis <- ssModel$R
lsel <- lower.tri(Rmis$values, TRUE)
mval <- Rmis$values[lsel]
mval[2] <- .05
mfre <- Rmis$free[lsel]
mfre[2] <- TRUE
mlab <- Rmis$labels[lsel]
Rmis <- mxMatrix("Symm", nrow(Rmis), ncol(Rmis), values=mval, labels=mlab, free=mfre, name=Rmis$name)
ssMis <- mxModel(ssModel, Rmis)

set.seed(101)
ssMisData <- mxGenerateData(ssMis, 200)

ssMisRun <- mxRun(mxModel(ssModel, mxData(ssMisData, 'raw')))

mi.mis <- mxMI(ssMisRun, full=FALSE)

if( any(mi.mis$MI > qchisq(p=1-0.01, df=1), na.rm=TRUE)){
	print("Large modification index found")
	# Grab the plus.param model that has the highest modification index
	newModelAttempt <- mi.mis$plusOneParamModels[[which.max(mi.mis$MI)]]
	newRunAttempt <- mxRun(newModelAttempt)
	mi.cor <- mxMI(newRunAttempt, full=FALSE)
	foundFirst <- TRUE
	if(  any(mi.cor$MI > qchisq(p=1-0.01, df=1), na.rm=TRUE ) ){
		print("Problem.  Model still may need modification")
		foundSecond <- TRUE
	} else{
		print("Success.  All modification indices are sufficiently small now.")
		foundSecond <- FALSE
	}
} else{
	foundFirst <- FALSE
}

omxCheckTrue(all( c(foundFirst, foundSecond) == c(TRUE, FALSE) ) )

prevMIsForSS <- c(0.000002, 10.827511, 1.664330, 2.062372, 1.211038,
	1.428606, 3.021569, 1.749314, 0.752319, 2.663316, 1.167776,
	0.032929, -0.186243)

miSSCheck <- mi.mis$MI[!is.na(mi.mis$MI)]
miSSCheck <- miSSCheck[-length(miSSCheck)]
prevMIsForSS <- prevMIsForSS[-length(prevMIsForSS)]

omxCheckCloseEnough(miSSCheck, prevMIsForSS, 0.01)

#------------------------------------------------------------------------------
# State space model check
# Similar to the above, but tests the respect of fixed parameter labels.
# When a fixed parameter occurs in more than one spot (and thus necessarily
#  has a label), it is freed by its label instead of by the matrix element.


# Create data based on state space model.
require(OpenMx)
nvar <- 5
varnames <- paste("x", 1:nvar, sep="")
ssModel <- mxModel(model="State Space Manual Example",
    mxMatrix("Full", 1, 1, TRUE, .3, name="A"),
    mxMatrix("Zero", 1, 1, name="B"),
    mxMatrix("Full", nvar, 1, TRUE, .6, name="C", dimnames=list(varnames, "F1")),
    mxMatrix("Zero", nvar, 1, name="D"),
    mxMatrix("Diag", 1, 1, FALSE, 1, name="Q"),
    mxMatrix("Diag", nvar, nvar, c(F, F, T, T, T), .2, name="R", labels=paste('resid', c(1, 1, 3:nvar), sep='')),
    mxMatrix("Zero", 1, 1, name="x0"),
    mxMatrix("Diag", 1, 1, FALSE, 1, name="P0"),
    mxMatrix("Zero", 1, 1, name="u"),
    mxExpectationStateSpace("A", "B", "C", "D", "Q", "R", "x0", "P0", "u"),
    mxFitFunctionML()
)


Rmis <- ssModel$R
lsel <- lower.tri(Rmis$values, TRUE)
mval <- Rmis$values[lsel]
mval[2] <- .05
mfre <- Rmis$free[lsel]
mfre[2] <- TRUE
mlab <- Rmis$labels[lsel]
Rmis <- mxMatrix("Symm", nrow(Rmis), ncol(Rmis), values=mval, labels=mlab, free=mfre, name=Rmis$name)
ssMis <- mxModel(ssModel, Rmis)

set.seed(101)
ssMisData <- mxGenerateData(ssMis, 200)

ssMisRun <- mxRun(mxModel(ssModel, mxData(ssMisData, 'raw')))

mi.mis <- mxMI(ssMisRun, full=FALSE)

if( any(mi.mis$MI > qchisq(p=1-0.01, df=1), na.rm=TRUE)){
	print("Large modification index found")
	# Grab the plus.param model that has the highest modification index
	newModelAttempt <- mi.mis$plusOneParamModels[[which.max(mi.mis$MI)]]
	newRunAttempt <- mxRun(newModelAttempt)
	mi.cor <- mxMI(newRunAttempt, full=FALSE)
	foundFirst <- TRUE
	if(  any(mi.cor$MI > qchisq(p=1-0.01, df=1), na.rm=TRUE ) ){
		print("Problem.  Model still may need modification")
		foundSecond <- TRUE
	} else{
		print("Success.  All modification indices are sufficiently small now.")
		foundSecond <- FALSE
	}
} else{
	foundFirst <- FALSE
}

omxCheckTrue(all( c(foundFirst, foundSecond) == c(TRUE, FALSE) ) )


miSSCheck <- mi.mis$MI[!is.na(mi.mis$MI)]
miSSCheck <- miSSCheck['resid1']
prevMIsForSS <- 3.909195

omxCheckCloseEnough(miSSCheck, prevMIsForSS, 0.01)



