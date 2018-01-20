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
# Program: MultiRegRaw-OpenMx100214.R
#  Author: Steven M. Boker
#    Date: Sun Feb 14 13:34:33 EST 2010
#
# This program fits a FIML multiple regression model to the 
#    multiData simulated data.
#
#
# ---------------------------------------------------------------------
# Revision History
#    -- Sun Feb 14 13:34:36 EST 2010
#      Created MultiRegRaw-OpenMx100214.R.
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

multiRegModel <- mxModel("FIML Multiple Regression of y on x1, x2, x3, and x4",
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
    mxPath(from="one", to=manifests, 
           arrows=1, free=TRUE, values=.1, 
           labels=c("MeanX1", "MeanX2", "MeanX3", "MeanX4", "MeanY")),
    mxData(observed=multiData1, type="raw")
    )

multiRegModelOut <- mxRun(multiRegModel, suppressWarnings=TRUE)

summary(multiRegModelOut)


# ----------------------------------
# check for correct values

expectVal <- c(0.04427, 0.306983, 0.398635, 0.471855, 1.136428, 
0.581104, 1.555603, 0.634979, 0.564913, 2.102767, 0.590538, 0.434355, 
0.650669, 2.552983, 0.533376, 0.984891, 1.974177, 3.073288, 4.055537, 
0.060049)

expectSE <- c(0.036967, 0.029548, 0.025327, 0.022115, 0.07186, 0.06488, 0.09838, 
0.074715, 0.084718, 0.132974, 0.080594, 0.091196, 0.107596, 0.161443, 
0.033734, 0.047675, 0.055779, 0.06485, 0.071456, 0.106581)

expectMin <- 7615.126

omxCheckCloseEnough(expectMin, multiRegModelOut$output$minimum, 0.001)

omxCheckCloseEnough(expectVal, multiRegModelOut$output$estimate, 0.002)

omxCheckCloseEnough(expectSE, 
    as.vector(multiRegModelOut$output[['standardErrors']]), 0.001)
