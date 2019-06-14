#
#   Copyright 2007-2019 by the individuals mentioned in the source code history
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

obsAcov <- structure(c(91.056, -23.916, -2.33, -27.78, -3.783, 4.485, 3.81,  -9.785, 
                       -23.916, 114.052, -34.385, -31.261, -4.954, 2.362, -3.783,  3.717, 
                       -2.33, -34.385, 68.399, 0.026, -15.058, -11.52, 8.942,  -7.985, -27.78,
                       -31.261, 0.026, 117.372, -38.064, -8.443, 3.329,  9.487, -3.783, -4.954,
                       -15.058, -38.064, 76.331, 12.763, -12.746,  -1.576, 4.485, 2.362, -11.52,
                       -8.443, 12.763, 121.513, -49.804,  -25.08, 3.81, -3.783, 8.942, 3.329,
                       -12.746, -49.804, 134.231,  -23.096, -9.785, 3.717, -7.985, 9.487, -1.576,
                       -25.08, -23.096,  181.501), .Dim = c(8L, 8L),
                     .Dimnames = list(c("xt1", "yt1",  "yt2", "zt1", "zt2", "poly_y_x", "poly_z_x", "poly_z_y"),
                                      c("xt1",  "yt1", "yt2", "zt1", "zt2", "poly_y_x", "poly_z_x", "poly_z_y" )))

obsThr <- matrix(c(
 0.02506887, 0.02506891, 0.1509692,
         NA, 0.84162123, 1.1749868),
	2, 3, byrow=TRUE)
dimnames(obsThr) <- list(NULL, letters[(27-nvar):26])

obsMns <- as.numeric(matrix(0, 1, 3))
names(obsMns) <- dimnames(obsCor)[[2]]

obsStats <- list(means=obsMns, cov=obsCor, thresholds=obsThr, acov=obsAcov, fullWeight=obsAcov)

obsWDat <- mxData(observed=rawData, type='raw', observedStats=obsStats)

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

omxCheckCloseEnough(mxEval(satCov, model=satwls2Run), obsStats$cov, 1e-4)
omxCheckCloseEnough(mxEval(thresholdMatrix, model=satwls2Run)[-2], obsStats$thresholds[-2], 1e-2)
omxCheckCloseEnough(mxEval(fitfunction, model=satwls2Run), 0, 1e-3)



#e <- c(1, .8, .8, 1, .8, 1., .1, .1, .8, .1, .8)
#o <- c(vech(testOld$observed), testOld$thresholds[-2])

#t(o-e) %*% testOld$acov %*% (o-e)  # Correct fit function value

#w <- o-e
#w[3] <- 0
#t(w) %*% testOld$acov %*% w  # Computed fit function value


