#
#   Copyright 2007-2021 by the individuals mentioned in the source code history
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
# Date: 2016.12.24
# Filename: jointFactorWls.R
# Purpose: Test joint ordinal and continuous weighted least squares.
#  In particular, I want to verify that for correctly specified models
#  full WLS gives unbiased estimates like diagonally weighted least
#  squares and unweighted least squares do.
# Note: This test is based on
#  inst/models/passing/jointFactorModelsTest.R
#------------------------------------------------------------------------------



#------------------------------------------------------------------------------
# Load package define a useful function


library(OpenMx)
library(testthat)

if (mxOption(NULL, "Default optimizer") == 'NPSOL') stop('SKIP')

rms <- function(x, y=NA){
	if(is.matrix(x) && is.vector(y) && nrow(x) == length(y)){
		sqrt(colMeans((x-y)^2))
	} else if(is.matrix(x) && length(y) == 1 && is.na(y)){
		apply(x, 2, FUN=rms, x=x)
	} else{
		sqrt(mean((x-y)^2))
	}
}



#------------------------------------------------------------------------------
# Get and Structure Data

data("jointdata", package ="OpenMx", verbose= TRUE)
jointData <- jointdata

# specify ordinal columns as ordered factors
jointData[,c(2,4,5)] <- mxFactor(jointData[,c(2,4,5)],
	levels=list(c(0,1), c(0, 1, 2, 3), c(0, 1, 2)))

if(1) {
  # Check joint WLS data generation
  set.seed(23)
  simData <- mxGenerateData(jointData)
  omxCheckTrue(dim(simData) == c(250, 5))
  omxCheckCloseEnough(cor(simData$z1, simData$z3), cor(jointData$z1, jointData$z3), 0.01)

  tabo <- table(simData$z2, simData$z4)
  tabe <- table(jointData$z2, jointData$z4)
  # Chi-square-ish test
  # That is, the vectorized joint distributions are near
  #  the expected value of the chi-square
  omxCheckTrue(sum((tabo-tabe)^2/tabe)/(4*2-1) < 1.6)

  omxCheckTrue(all.equal(sapply(jointData, levels), sapply(simData, levels)))
}


#------------------------------------------------------------------------------
# Model pre-setting

satCov <- mxMatrix("Symm", 5, 5,
	free=TRUE, values=diag(5), name="C")
satCov$free[2,2] <- FALSE
satCov$free[4,4] <- FALSE
satCov$free[5,5] <- FALSE

loadings <- mxMatrix("Full", 1, 5,
	free=TRUE, values=1, name="L", lbound=0, labels=paste0('l',1:5))
loadings$ubound[1,5] <- 2

resid <- mxMatrix("Diag", 5, 5,
	free=c(TRUE, FALSE, TRUE, FALSE, FALSE), values=.5, name="U")

means <- mxMatrix("Full", 1, 5,
	free=c(TRUE, FALSE, TRUE, FALSE, FALSE), values=0, name="M")

thresh <- mxMatrix("Full", 3, 3, FALSE, 0, name="T")

thresh$free[,1] <- c(TRUE, FALSE, FALSE)
thresh$values[,1] <- c(0, NA, NA)
thresh$labels[,1] <- c("z2t1", NA, NA)

thresh$free[,2] <- TRUE
thresh$values[,2] <- c(-1, 0, 1)
thresh$labels[,2] <- c("z4t1", "z4t2", "z4t3")

thresh$free[,3] <- c(TRUE, TRUE, FALSE)
thresh$values[,3] <- c(-1, 1, NA)
thresh$labels[,3] <- c("z5t1", "z5t2", NA)
colnames(thresh) <- paste0('z', c(2,4,5))

#------------------------------------------------------------------------------
# Model definition

# ML form
jointModel1 <- mxModel("ContinuousOrdinalData",
				mxData(jointData, "raw"),
				loadings, resid, means, thresh,
			mxAlgebra(t(L) %*% L + U, name="C"),
			mxFitFunctionML(),
			mxExpectationNormal("C", "M",
				dimnames=names(jointData),
				thresholds="T",
				threshnames=c("z2", "z4", "z5"))
			)

ramModel1 <- jointRAM <- mxModel(
  "JointRAM", type="RAM", thresh,
  manifestVars = paste0('z', 1:5),
  latentVars = 'G',
  mxData(jointData, "raw"),
  mxPath('one', paste0('z', c(1,3)), free=TRUE),
  mxPath(paste0('z', c(1,3)), arrows=2, free=TRUE, values=.5),
  mxPath(paste0('z', c(2,4,5)), arrows=2, free=FALSE, values=.5),
  mxPath('G', arrows=2, values=1, free=FALSE),
  mxPath('G', paste0('z', 1:5),
         free=TRUE, values=1, lbound=0))

ramModel1$expectation$thresholds <- 'T'

# Create WLS Data
wd <- mxData(jointData, "raw")
dd <- mxData(jointData, "raw")
ud <- mxData(jointData, "raw")

# WLS form(s) of model
jointWlsModel <- mxModel(jointModel1, name='wlsModel', wd, mxFitFunctionWLS('WLS'))
jointDlsModel <- mxModel(jointModel1, name='dlsModel', dd, mxFitFunctionWLS('DWLS'))
jointUlsModel <- mxModel(jointModel1, name='ulsModel', ud, mxFitFunctionWLS('ULS'))

ramWlsModel <- mxModel(ramModel1, name='wlsModel', mxFitFunctionWLS())

# Run 'em
jointWlsResults <- mxRun(jointWlsModel)
jointDlsResults <- mxRun(jointDlsModel)
jointUlsResults <- mxRun(jointUlsModel)

expect_equal(diag(jointWlsResults$data$observedStats$asymCov),
             diag(jointDlsResults$data$observedStats$asymCov))

ramWlsResults <- mxRun(ramWlsModel)

jointResults1 <- mxRun(jointModel1, suppressWarnings = TRUE)
summary(jointResults1)

ramResult1 <- mxRun(ramModel1)
summary(ramResult1)

omxCheckCloseEnough(coef(jointWlsResults) - coef(ramWlsResults),
                    rep(0,15), 1e-3)

expect_equal(jointDlsResults$output$fit, 5.061, .1)
expect_equal(summary(jointDlsResults)$Chi, 10.126, .1)
expect_equal(summary(jointDlsResults)$ChiDoF, 5)

noOptResult <- mxRun(jointDlsResults, useOptimizer = FALSE)
expect_true(is.null(summary(noOptResult)$Chi))

expect_equal(length(mxGetExpected(jointWlsResults, 'standvector')),
             summary(jointWlsResults)$observedStatistics)

#------------------------------------------------------------------------------
# Compare ML and WLS estimates

round(cmp <- cbind(ML=coef(jointResults1), WLS=coef(jointWlsResults), DLS=coef(jointDlsResults), ULS=coef(jointUlsResults)), 3)

plot(cmp[1:5,1], cmp[1:5,2])
abline(0, 1)

print(rms(cmp))
omxCheckCloseEnough(vechs(rms(cmp)), c(0.01, 0.02, 0.029, 0.027, 0.036, 0.01), .002)

round(seCmp <- cbind(ML=jointResults1$output$standardErrors,
                     WLS=jointWlsResults$output$standardErrors,
                     DLS=jointDlsResults$output$standardErrors,
                     ULS=jointUlsResults$output$standardErrors), 3)

omxCheckCloseEnough(cor(seCmp)[1,2:4], rep(1,3), .04)

se1 <- c(0.065, 0.113, 0.072, 0.106, 0.068, 0.064, 0.067,  0.061,
         0.061, 0.074, 0.083, 0.075, 0.098, 0.067, 0.06)
omxCheckCloseEnough(c(jointDlsResults$output$standardErrors), se1, .01)

#------------------------------------------------------------------------------
# Create and compare saturated models

jointModel2 <- mxModel("ContinuousOrdinalData",
				mxData(jointData, "raw"),
				satCov, means, thresh,
				mxFitFunctionML(),
				mxExpectationNormal("C", "M",
					dimnames=names(jointData),
					thresholds="T",
					threshnames=c("z2", "z4", "z5"))
				)

jointResults2 <- mxRun(jointModel2, suppressWarnings = TRUE)
summary(jointResults2)


ref <- mxRefModels(jointModel1, run=TRUE)

(sref <- summary(jointResults1, refModels=ref))

(shan <- summary(jointResults1, SaturatedLikelihood=-2*logLik(jointResults2), SaturatedDoF=summary(jointResults2)$degreesOfFreedom))

(swls <- summary(jointWlsResults))

# Compare chi-squared degrees of freedom
omxCheckTrue(sref$ChiDoF == shan$ChiDoF)
omxCheckTrue(shan$ChiDoF == swls$ChiDoF)


# Compare Chi-squared values too
omxCheckCloseEnough(sref$Chi, as.numeric(shan$Chi), 1e-3)
omxCheckWithinPercentError(as.numeric(shan$Chi), swls$Chi, 28)



#------------------------------------------------------------------------------
# Check that ML saturated model estimates are close
#  to the WLS saturated model estimates.

obsStats <- jointWlsResults$data$observedStats

expect_error(mxGetExpected(jointResults2, 'covariance', subname="xyz"),
             "not found")

mxGetExpected(jointResults2, 'covariance')
obsStats$cov

mxGetExpected(jointResults2, 'means')
obsStats$means

mxGetExpected(jointResults2, 'thresholds')
obsStats$thresholds

ml.sat <- mxGetExpected(jointResults2, 'vector')
wls.sat <- c(vech(obsStats$cov), obsStats$means, na.omit(c(obsStats$thresholds)))

omxCheckTrue(rms(ml.sat, wls.sat) < .01)
omxCheckCloseEnough(ml.sat, wls.sat, .03) #could adjust to 0.009


#------------------------------------------------------------------------------
# Check that observed statistics are correct

numCovs <- 5*4/2
numMeans <- 2
numVaris <- 2
numLevel <- sapply(sapply(jointData, levels), length)
numThres <- sum((numLevel - 1) * (numLevel > 0))
numObsSummaryStats <- numCovs + numMeans + numVaris + numThres

omxCheckEquals(summary(jointWlsResults)$observedStatistics, numObsSummaryStats)
omxCheckEquals(summary(jointDlsResults)$observedStatistics, numObsSummaryStats)
omxCheckEquals(summary(jointUlsResults)$observedStatistics, numObsSummaryStats)

omxCheckEquals(summary(jointResults1)$observedStatistics, prod(dim(jointData)))


#------------------------------------------------------------------------------
# Regularization Test

regTest <- mxModel(jointWlsModel,
        mxPenaltyLASSO(paste0('l', 1:5), name="lasso"),
        mxMatrix(nrow=1, ncol=1, free=TRUE, values=1, labels="lambda"))
regTest <- mxRun(regTest)
expect_equivalent(regTest$fitfunction$result[1,1],
                   mxEval(fitfunction, regTest)[1,1])
expect_equivalent(regTest$lasso$result[1,1],
             mxEval(lasso, regTest))
expect_equivalent(regTest$lasso$result[1,1],
             mxEval(lasso, regTest, compute=TRUE))
