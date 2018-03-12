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
# Program: UnivariateRaw-OpenMx100214.R
#  Author: Steven M. Boker
#    Date: Sun Feb 14 13:18:47 EST 2010
#
# This program fits a FIML univariate model to the 
#     multiData simulated data.
#
#
# ---------------------------------------------------------------------
# Revision History
#    -- Sun Feb 14 13:18:51 EST 2010
#      Created UnivariateRaw-OpenMx100214.R.
#
# ---------------------------------------------------------------------

# ----------------------------------
# Read libraries and set options.

require(OpenMx)

# ----------------------------------
# Read the data and print descriptive statistics.

data(multiData1)

# ----------------------------------
# Build an OpenMx univariate regression model using y and x1

manifests <- c("x1", "y")

uniRegModelRaw <- mxModel("ur",
    type="RAM",
    manifestVars=manifests,
    mxPath(from="x1", to="y", arrows=1, 
           free=TRUE, values=.2, labels="b1"),
    mxPath(from=manifests, 
           arrows=2, free=TRUE, values=.8, 
           labels=c("VarX1", "VarE")),
    mxPath(from="one", to=manifests, 
           arrows=1, free=TRUE, values=.1, 
           labels=c("MeanX1", "MeanY")),
    mxData(observed=multiData1, type="raw")
    )

uniRegModelRawOut <- mxRun(uniRegModelRaw)

summary(uniRegModelRawOut)


#---------------------
# check values: uniRegModelRawOut

expectVal <- c(0.669179, 1.13643, 1.647629, 0.984894, 3.189368)

expectSE <-c(0.053849, 0.071873, 0.104204, 0.047674, 0.078154)

expectMin <- 3151.492

omxCheckCloseEnough(expectVal, uniRegModelRawOut$output$estimate, 0.001)

omxCheckCloseEnough(expectSE, 
    as.vector(uniRegModelRawOut$output[['standardErrors']]), 0.001)

omxCheckCloseEnough(expectMin, uniRegModelRawOut$output$minimum, 0.001)

# -----------------

uniRegModelRaw$fitfunction$rowDiagnostics <- TRUE
urm1 <- mxRun(uniRegModelRaw)
lk <- attr(urm1$fitfunction$result, 'likelihoods')

omxCheckCloseEnough(urm1$output$fit, expectMin, 1e-3)
omxCheckCloseEnough(-2 * sum(log(lk)), expectMin, 1e-3)

algModel <- mxModel("wrap",
	mxModel(uniRegModelRaw, mxFitFunctionML(vector=TRUE, rowDiagnostics=TRUE)),
	mxAlgebra(-2*sum(log(ur.fitfunction)), 'm2ll'),
	mxFitFunctionAlgebra('m2ll'))

urm2 <- mxRun(algModel)

omxCheckCloseEnough(urm2$output$fit, expectMin, 1e-3)

omxCheckEquals(max(abs(urm2$ur$fitfunction$info$likelihoods -
                         urm2$ur$fitfunction$result)),0)
