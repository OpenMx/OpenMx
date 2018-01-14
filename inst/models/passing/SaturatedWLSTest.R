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

#------------------------------
# Author: Michael D. Hunter
# Date: 2014.03.27
# Filename: SaturatedWLSTest.R
# Purpose: Test the WLS fit function
#  in a saturated model.
#------------------------------


#------------------------------
# Load packages and old WLS data

require(mvtnorm)
require(OpenMx)

#load("wlsTest1.RData")

f1 <- ordered(c("a","b"), levels=c("a","b"))
f1 <- ordered(f1, levels=c("b"), exclude="a")

f2 <- mxFactor(c("a","b"), levels=c("a","b"))
f2 <- omxCheckWarning(mxFactor(f2, levels=c("b"), exclude="a"), NA)

omxCheckEquals(f1[2], f2[2])
omxCheckEquals(is.na(f1), is.na(f2))
omxCheckEquals(levels(f1), levels(f2))

#------------------------------
# Generate data

nvar <- 3
k    <- .5
sigma <- matrix(k, nvar, nvar)
diag(sigma) <- 1

set.seed(943)


cDat    <- rmvnorm(100, sigma=sigma)
rawData <- (cDat>0) + (cDat>1)
rawData[rawData[,1]>1,1] <- 1

cCor <- cor(cDat)
oCor <- cor(rawData)

rawData <- data.frame(
	mxFactor(as.data.frame(rawData[,1]),      0:1),
	mxFactor(as.data.frame(rawData[,2:nvar]), 0:2)
	)
names(rawData) <- letters[(27-nvar):26]

cDat <- data.frame(cDat)

#------------------------------
obsCor <- matrix(c(
	1.0000000, 0.6032634, 0.5532893,
	0.6032634, 1.0000000, 0.5160187,
	0.5532893, 0.5160187, 1.0000000),
	3, 3, byrow=TRUE)
dimnames(obsCor) <- list(letters[(27-nvar):26], letters[(27-nvar):26])
obsCor <- round(obsCor, 6)

obsAcov <- matrix(c(
    0,   0.000000,   0.000000,    0,   0.000000,    0,   0.000000,   0.000000,   0.0000000,   0.0000000,   0.000000,
    0, 124.645281, -44.386060,    0, -13.665562,    0,  -3.486089,   7.949569,   9.9267071,  -5.8507369,  -8.066406,
    0, -44.386060,  97.021057,    0,  -9.126442,    0,  12.745952,  -8.916994,  -9.9068741,   3.5362051,  -7.396720,
    0,   0.000000,   0.000000,    0,   0.000000,    0,   0.000000,   0.000000,   0.0000000,   0.0000000,   0.000000,
    0, -13.665562,  -9.126442,    0, 116.569069,    0, -19.654759,  -2.427452,  -2.0451332,   7.8230997,  11.144028,
    0,   0.000000,   0.000000,    0,   0.000000,    0,   0.000000,   0.000000,   0.0000000,   0.0000000,   0.000000,
    0,  -3.486089,  12.745952,    0, -19.654759,    0,  92.629655, -19.657405, -14.8964215, -29.3524731,   1.615608,
    0,   7.949569,  -8.916994,    0,  -2.427452,    0, -19.657405,  98.441444, -31.0103427, -10.8982547,  -9.676036,
    0,   9.926707,  -9.906874,    0,  -2.045133,    0, -14.896422, -31.010343,  71.4383048,  -0.9634742,  -7.085017,
    0,  -5.850737,   3.536205,    0,   7.823100,    0, -29.352473, -10.898255,  -0.9634742,  91.5291529, -21.750644,
    0,  -8.066406,  -7.396720,    0,  11.144028,    0,   1.615608,  -9.676036,  -7.0850173, -21.7506439,  51.090309
	),
	11, 11, byrow=TRUE)
obsAcov <- round(obsAcov, 6)
obsAcov <- cbind(obsAcov[,1:6], matrix(0, 11, 3), obsAcov[,7:11])
obsAcov <- rbind(obsAcov[1:6,], matrix(0, 3, 14), obsAcov[7:11,])

obsThr <- matrix(c(
 0.02506887, 0.02506891, 0.1509692,
         NA, 0.84162123, 1.1749868),
	2, 3, byrow=TRUE)
dimnames(obsThr) <- list(NULL, letters[(27-nvar):26])

obsMns <- as.numeric(matrix(0, 1, 3))
names(obsMns) <- dimnames(obsCor)[[2]]

obsWDat <- mxData(observed=obsCor, type='acov', means=obsMns, thresholds=obsThr, acov=obsAcov, fullWeight=obsAcov, numObs=nrow(cDat))

#------------------------------
# Make WLS saturated model

theDims <- list(c('x', 'y', 'z'), c('x', 'y', 'z'))

satwls2 <- mxModel(name="ExpNormWLSSat",
	mxMatrix("Symm", 3, 3, values=c(1, .8, .8, 1, .8, 1), free=c(FALSE,TRUE,TRUE,FALSE,TRUE,FALSE), dimnames=theDims, name="satCov"),
	mxMatrix("Full", 2, 3, values=c(0, NA, 0, .8, 0, .8), free=c(TRUE,FALSE,TRUE,TRUE,TRUE,TRUE), name="thresholdMatrix"),
	mxMatrix('Zero', nrow=1, ncol=3, name='meansMatrix'),
	obsWDat,
	mxFitFunctionWLS(),
	mxExpectationNormal(covariance="satCov", means='meansMatrix', thresholds="thresholdMatrix", dimnames=theDims[[1]])
)

#satwls2 <- mxOption(satwls2, "Calculate Hessian", "No")
#satwls2 <- mxOption(satwls2, "Standard Errors", "No")
#satwls2 <- mxOption(satwls2, "Major iterations", 1)






satwls2Run <- mxRun(satwls2)


#------------------------------
# Compare saturated model estimates to acov data

omxCheckCloseEnough(mxEval(satCov, model=satwls2Run), obsWDat$observed, 1e-4)
omxCheckCloseEnough(mxEval(thresholdMatrix, model=satwls2Run)[-2], obsWDat$thresholds[-2], 1e-2)
omxCheckCloseEnough(mxEval(fitfunction, model=satwls2Run), 0, 1e-3)



#e <- c(1, .8, .8, 1, .8, 1., .1, .1, .8, .1, .8)
#o <- c(vech(testOld$observed), testOld$thresholds[-2])

#t(o-e) %*% testOld$acov %*% (o-e)  # Correct fit function value

#w <- o-e
#w[3] <- 0
#t(w) %*% testOld$acov %*% w  # Computed fit function value


