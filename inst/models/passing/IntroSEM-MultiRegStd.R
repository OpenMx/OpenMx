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
# Program: MultiRegStd-OpenMx100214.R
#  Author: Steven M. Boker
#    Date: Sun Feb 14 12:25:16 EST 2010
#
# This program fits a multiple regression model to the 
#    multiData simulated data.
#
#
# ---------------------------------------------------------------------
# Revision History
#    -- Sun Feb 14 12:25:21 EST 2010
#      Created MultiRegStd-OpenMx100214.R.
#
# ---------------------------------------------------------------------

# ----------------------------------
# Read libraries and set options.

require(OpenMx)

# ----------------------------------
# Read the data and print descriptive statistics.

data(multiData1)

# ----------------------------------
# Build an OpenMx multiple regression model using y and x1

predictors <- c("x1", "x2", "x3", "x4")
outcomes <- c("y")
manifests <- names(multiData1)
multiData1Cov <- cov(multiData1)

multiRegModel <- mxModel("Multiple Regression of y on x1, x2, x3, and x4",
    type="RAM",
    manifestVars=manifests,
    mxPath(from=predictors, to=outcomes, 
           arrows=1, 
           free=TRUE, values=.2, 
           labels=c("b1", "b2", "b3", "b4")),
    mxPath(from=outcomes, 
           arrows=2, 
           free=TRUE, values=.8, 
           labels=c("VarE")),
    mxPath(from=predictors, to=predictors,
#           arrows=2, all=TRUE,
           arrows=2, connect="unique.pairs",
           free=TRUE, values=.2),
    mxPath(from=manifests, 
           arrows=2, 
           free=TRUE, values=.8, 
           labels=c("VarX1", "VarX2", "VarX3", "VarX4", "VarE")),
    mxData(observed=multiData1Cov, type="cov", numObs=500)
    )

multiRegModelOut <- mxRun(multiRegModel, suppressWarnings=TRUE)

summary(multiRegModelOut)


# ----------------------------------
# check for correct values

expectVal <- c(0.04427, 0.30698, 0.39864, 0.47186, 1.13643, 0.58111,
  1.5556, 0.63498, 0.56491, 2.10277, 0.59054, 0.43436, 0.65067,  2.55298, 0.53338)

expectSE <- c(0.037002, 0.029577, 0.025352, 0.022138, 0.072093, 0.065091, 
0.098683, 0.074968, 0.084998, 0.133396, 0.080868, 0.091494, 0.107955, 
0.161951, 0.033835)

# cat(deparse(round(multiRegModelOut$output$estimate, 5)))
omxCheckCloseEnough(expectVal, multiRegModelOut$output$estimate, 0.001)

omxCheckCloseEnough(expectSE, 
    as.vector(multiRegModelOut$output[['standardErrors']]), 0.001)

omxCheckCloseEnough(3020.43369, multiRegModelOut$output$minimum, 0.001)


    

