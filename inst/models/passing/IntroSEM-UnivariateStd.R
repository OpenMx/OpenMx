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
# Program: UnivariateStd-OpenMx100214.R
#  Author: Steven M. Boker
#    Date: Sun Feb 14 12:13:20 EST 2010
#
# This program fits a univariate model to the multiData simulated data.
#
#
# ---------------------------------------------------------------------
# Revision History
#    -- Sun Feb 14 12:13:16 EST 2010
#      Created UnivariateStd-OpenMx100214.R.
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
multiData1Cov <- cov(multiData1[,c(1,5)])

uniRegModel <- mxModel("Univariate Regression of y on x1",
    type="RAM",
    manifestVars=manifests,
    mxPath(from="x1", to="y", arrows=1, 
           free=TRUE, values=.2, labels="b1"),
    mxPath(from=manifests, arrows=2, 
           free=TRUE, values=.8, labels=c("VarX1", "VarE")),
    mxData(observed=multiData1Cov, type="cov", numObs=500)
    )

uniRegModelOut <- mxRun(uniRegModel, suppressWarnings=TRUE)

summary(uniRegModelOut)

#---------------------
# check values: uniRegModelOut

expectVal <- c(0.66918, 1.13643, 1.64763)

expectSE <-c(0.053902, 0.07209, 0.104518)

# cat(deparse(round(uniRegModelOut$output$estimate, 5)))
omxCheckCloseEnough(expectVal, uniRegModelOut$output$estimate, 0.001)

omxCheckCloseEnough(expectSE, 
    as.vector(uniRegModelOut$output[['standardErrors']]), 0.001)

omxCheckCloseEnough(1313.6145, uniRegModelOut$output$minimum, 0.001)

omxCheckEquals(uniRegModelOut$output$fitUnits, "-2lnL")

