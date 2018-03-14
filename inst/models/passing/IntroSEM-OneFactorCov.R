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


# ---------------------------------------------------------------------
# Program: OneFactorCov-OpenMx100221.R
#  Author: Steven M. Boker
#    Date: Sun Feb 21 13:39:29 EST 2010
#
# This program fits a covariance single factor model to the 
#     factorExample1.csv simulated data.
#
#
# ---------------------------------------------------------------------
# Revision History
#    -- Sun Feb 21 13:39:32 EST 2010
#      Created OneFactorCov-OpenMx100221.R.
#
# ---------------------------------------------------------------------

# ----------------------------------
# Read libraries and set options.

require(OpenMx)

# ----------------------------------
# Read the data and print descriptive statistics.

data(factorExample1)

# ----------------------------------
# Build an OpenMx single factor covariance model with fixed variance

indicators <- names(factorExample1)
latents <- c("F1")
loadingLabels <- paste("b_", indicators, sep="")
uniqueLabels <- paste("U_", indicators, sep="")
meanLabels <- paste("M_", indicators, sep="")
factorVarLabels <- paste("Var_", latents, sep="")

oneFactorCov1 <- mxModel("Single Factor Covariance Model with Fixed Variance",
    type="RAM",
    manifestVars=indicators,
    latentVars=latents,
    mxPath(from=latents, to=indicators, 
#           arrows=1, all=TRUE, 
           arrows=1, connect="all.pairs", 
           free=TRUE, values=.2, 
           labels=loadingLabels),
    mxPath(from=indicators, 
           arrows=2, 
           free=TRUE, values=.8, 
           labels=uniqueLabels),
    mxPath(from=latents,
           arrows=2, 
           free=FALSE, values=1, 
           labels=factorVarLabels),
    mxData(observed=cov(factorExample1), type="cov", numObs=500)
    )

oneFactorCov1Out <- mxRun(oneFactorCov1, suppressWarnings=TRUE)

summary(oneFactorCov1Out)

omxCheckError(mxBootstrap(oneFactorCov1Out),
              "MxComputeBootstrap: data 'Single Factor Covariance Model with Fixed Variance.data' of type 'cov' cannot have row weights")

# ----------------------------------
# check for correct values

expectVal <- c(0.68396, 0.32482, 0.10887, 0.47441, 0.60181, 1.12064,
1.25934, 0.64739, 0.71873, 0.3528, 0.17619, 0.19354, 0.79988, 0.63306,
0.36763, 0.34024, 0.23404, 0.85441)

expectSE <- c(0.035244, 0.022431, 0.020808, 0.044662, 0.042295, 0.045789, 
0.048931, 0.03064, 0.049368, 0.02492, 0.01197, 0.01234, 0.052165, 
0.042853, 0.032175, 0.034939, 0.017353, 0.057948)

# cat(deparse(round(oneFactorCov1Out$output$estimate, 5)))
omxCheckCloseEnough(expectVal, oneFactorCov1Out$output$estimate, 0.001)

omxCheckCloseEnough(expectSE, 
    as.vector(oneFactorCov1Out$output[['standardErrors']]), 0.001)

omxCheckCloseEnough(1435.94, oneFactorCov1Out$output$minimum, 0.001)


# ----------------------------------
# Check that WLS can be swapped into a path model just by
#  adding a WLS data set.
oneFactorCovWLS <- mxModel(oneFactorCov1Out, name='WLS',
	mxDataWLS(factorExample1)
)

oneFactorCovWLSOut <- mxRun(oneFactorCovWLS)

# WLS estimates are close to ML estimates
rms <- function(x, y){sqrt(mean((x-y)^2))}
omxCheckTrue(rms(expectVal, omxGetParameters(oneFactorCovWLSOut)) < .25)
omxCheckTrue(rms(expectSE, summary(oneFactorCovWLSOut)$parameters[,6]) < .025)


# Swap the data back to cov and make sure the fit function automatically adjusts
oneFactorCovML <- mxModel(oneFactorCovWLSOut, name='ML',
	mxData(observed=cov(factorExample1), type="cov", numObs=500)
)

oneFactorCovMLOut <- mxRun(oneFactorCovML, suppressWarnings=TRUE)
omxCheckCloseEnough(expectVal, oneFactorCovMLOut$output$estimate, 0.001)


# ----------------------------------
# Build an OpenMx single factor covariance model with fixed loading


oneFactorCov2 <- mxModel("Single Factor Covariance Model with Fixed Loading",
    type="RAM",
    manifestVars=indicators,
    latentVars=latents,
    mxPath(from=latents, to=indicators, 
           arrows=1, connect="all.pairs",
           free=TRUE, values=.2, 
           labels=loadingLabels),
    mxPath(from=indicators, 
           arrows=2, 
           free=TRUE, values=.8, 
           labels=uniqueLabels),
    mxPath(from=latents,
           arrows=2,
           free=TRUE, values=1, 
           labels=factorVarLabels),
    mxPath(from=latents, to=c("x1"),
           arrows=1, 
           free=FALSE, values=1),
    mxData(observed=cov(factorExample1), type="cov", numObs=500)
    )

oneFactorCov2Out <- mxRun(oneFactorCov2, suppressWarnings=TRUE)

summary(oneFactorCov2Out)

# ----------------------------------
# check for correct values

expectVal <- c(0.47491, 0.15917, 0.69363, 0.87989, 1.63847, 1.84125,
0.94654, 1.05084, 0.3528, 0.17619, 0.19354, 0.79988, 0.63306, 0.36763,
0.34024, 0.23404, 0.85441, 0.4678)

expectSE <- c(0.035182, 0.030686, 0.067917, 0.066228, 0.078835, 0.085033, 
0.051178, 0.077348, 0.02492, 0.01197, 0.01234, 0.052166, 0.042854, 
0.032175, 0.034939, 0.017353, 0.057948, 0.048261)

# cat(deparse(round(oneFactorCov2Out$output$estimate, 5)))
omxCheckCloseEnough(expectVal, oneFactorCov2Out$output$estimate, 0.001)

omxCheckCloseEnough(expectSE, 
    as.vector(oneFactorCov2Out$output[['standardErrors']]), 0.001)

omxCheckCloseEnough(1435.9409, oneFactorCov2Out$output$minimum, 0.001)


