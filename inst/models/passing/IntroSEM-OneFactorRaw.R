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
# Program: OneFactorRaw-OpenMx100221.R
#  Author: Steven M. Boker
#    Date: Sun Feb 21 13:23:40 EST 2010
#
# This program fits a FIML single factor model to the 
#     factorExample1.csv simulated data.
#
#
# ---------------------------------------------------------------------
# Revision History
#    -- Sun Feb 21 13:23:43 EST 2010
#      Created OneFactorRaw-OpenMx100221.R.
#
# ---------------------------------------------------------------------

# ----------------------------------
# Read libraries and set options.

require(OpenMx)

# ----------------------------------
# Read the data and print descriptive statistics.

data(factorExample1)

# ----------------------------------
# Build an OpenMx single factor FIML model with fixed variance

indicators <- names(factorExample1)
latents <- c("F1")
loadingLabels <- paste("b_", indicators, sep="")
uniqueLabels <- paste("U_", indicators, sep="")
meanLabels <- paste("M_", indicators, sep="")
factorVarLabels <- paste("Var_", latents, sep="")

oneFactorRaw1 <- mxModel("Single Factor FIML Model with Fixed Variance",
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
    mxPath(from="one", to=indicators, 
           arrows=1, free=TRUE, values=.1, 
           labels=meanLabels),
    mxData(observed=factorExample1, type="raw")
    )

oneFactorRaw1Out <- mxRun(oneFactorRaw1, suppressWarnings=TRUE)

summary(oneFactorRaw1Out)

# ----------------------------------
# check for correct values

expectVal <- c(0.683956, 0.32482, 0.108867, 0.474409, 0.601804, 
1.120639, 1.259331, 0.647393, 0.718727, 0.352796, 0.176193, 0.193536, 
0.799875, 0.633057, 0.367627, 0.340238, 0.234038, 0.854411, -0.011613, 
-0.006823, 0.023961, -0.031357, -0.045482, -0.091784, -0.067323, 
-0.03902, -0.059997)

expectSE <- c(0.035171, 0.022385, 0.020766, 0.044572, 0.042209, 0.045695, 
0.048829, 0.030576, 0.049267, 0.024845, 0.011934, 0.012303, 0.05201, 
0.042727, 0.032079, 0.034835, 0.017301, 0.057773, 0.04048, 0.023723, 
0.020266, 0.045256, 0.044588, 0.056927, 0.062004, 0.036112, 0.052335
)

expectMin <- 9706.388

omxCheckCloseEnough(expectVal, oneFactorRaw1Out$output$estimate, 0.001)

omxCheckCloseEnough(expectSE, 
    as.vector(oneFactorRaw1Out$output[['standardErrors']]), 0.001)

omxCheckCloseEnough(expectMin, oneFactorRaw1Out$output$minimum, 0.001)


# ----------------------------------
# Build an OpenMx single factor FIML model with fixed loading

indicators <- names(factorExample1)
latents <- c("F1")
loadingLabels <- paste("b_", indicators, sep="")
uniqueLabels <- paste("U_", indicators, sep="")
meanLabels <- paste("M_", indicators, sep="")
factorVarLabels <- paste("Var_", latents, sep="")

oneFactorRaw2 <- mxModel("Single Factor FIML Model with Fixed Loading",
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
           free=TRUE, values=1, 
           labels=factorVarLabels),
    mxPath(from=latents, to=c("x1"),
           arrows=1, 
           free=FALSE, values=1),
    mxPath(from="one", to=indicators, 
           arrows=1, free=TRUE, values=.1, 
           labels=meanLabels),
    mxData(observed=factorExample1, type="raw")
    )

oneFactorRaw2Out <- mxRun(oneFactorRaw2, suppressWarnings=TRUE)

summary(oneFactorRaw2Out)



# ----------------------------------
# check for correct values

expectVal <- c(0.474912, 0.159172, 0.693622, 0.879884, 1.638461, 
1.841241, 0.946538, 1.050835, 0.352795, 0.176193, 0.193536, 0.799875, 
0.633057, 0.367627, 0.340238, 0.234038, 0.854411, 0.467801, -0.011613, 
-0.006823, 0.023961, -0.031357, -0.045482, -0.091784, -0.067324, 
-0.039021, -0.059997)

expectSE <- c(0.035139, 0.030654, 0.067841, 0.066148, 0.078732, 0.084924, 
0.051113, 0.077259, 0.024844, 0.011934, 0.012303, 0.052011, 0.042727, 
0.03208, 0.034835, 0.017301, 0.057779, 0.048107, 0.040474, 0.023721, 
0.020265, 0.045254, 0.044584, 0.056914, 0.061993, 0.036106, 0.052327
)

expectMin <- 9706.388

omxCheckCloseEnough(expectVal, oneFactorRaw2Out$output$estimate, 0.001)

omxCheckCloseEnough(expectSE, 
    as.vector(oneFactorRaw2Out$output[['standardErrors']]), 0.001)

omxCheckCloseEnough(expectMin, oneFactorRaw2Out$output$minimum, 0.001)




