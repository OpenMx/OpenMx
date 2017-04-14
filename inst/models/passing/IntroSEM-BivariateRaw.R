#
#   Copyright 2007-2017 The OpenMx Project
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
# Program: BivariateRaw-OpenMx100214.R
#  Author: Steven M. Boker
#    Date: Sun Feb 14 13:27:16 EST 2010
#
# This program fits a FIML bivariate model to the 
#     multiData simulated data.
#
#
# ---------------------------------------------------------------------
# Revision History
#    -- Sun Feb 14 13:27:13 EST 2010
#      Created BivariateRaw-OpenMx100214.R.
#    -- Sun Aug 29 2010
#      Formatted for OpenMx Test Suite
# ---------------------------------------------------------------------

# ----------------------------------
# Read libraries and set options.

require(OpenMx)

# ----------------------------------
# Read the data and print descriptive statistics.

data(multiData1)

# ----------------------------------
# Build an OpenMx bivariate regression model using y and x1

manifests <- c("x1", "x2", "y")

biRegModelRaw <- mxModel("FIML Bivariate Regression of y on x1 and x2",
    type="RAM",
    manifestVars=manifests,
    mxPath(from=c("x1","x2"), to="y", 
           arrows=1, 
           free=TRUE, values=.2, labels=c("b1", "b2")),
    mxPath(from=manifests, 
           arrows=2, 
           free=TRUE, values=.8, 
           labels=c("VarX1", "VarX2", "VarE")),
    mxPath(from="x1", to="x2",
           arrows=2, 
           free=TRUE, values=.2, 
           labels=c("CovX1X2")),
    mxData(observed=multiData1, type="raw")
    )

omxCheckError(mxRun(biRegModelRaw, suppressWarnings=TRUE),
"The job for model 'FIML Bivariate Regression of y on x1 and x2' exited abnormally with the error message: FIML Bivariate Regression of y on x1 and x2.expectation: raw data observed but no expected means vector was provided. Add something like mxPath(from = 'one', to = manifests) to your model.")

biRegModelRaw <- mxModel(
  biRegModelRaw,
  mxPath(from="one", to=manifests, 
         arrows=1, free=TRUE, values=.1, 
         labels=c("MeanX1", "MeanX2", "MeanY")))

biRegModelRawOut <- mxRun(biRegModelRaw)

summary(biRegModelRawOut)


# ----------------------------------
# check for correct values

expectVal <- c(0.4479, 0.4328, 1.1364, 0.5811, 1.5556, 1.412, 0.9849, 
1.9741, 2.5529)

expectSE <- c(0.0554, 0.0474, 0.0719, 0.0649, 0.0984, 0.0893, 0.0477, 0.0558, 
0.1004)

expectMin <- 4608.207

omxCheckCloseEnough(expectVal, biRegModelRawOut$output$estimate, 0.001)

omxCheckCloseEnough(expectSE, 
    as.vector(biRegModelRawOut$output[['standardErrors']]), 0.001)

omxCheckCloseEnough(expectMin, biRegModelRawOut$output$minimum, 0.001)

omxCheckCloseEnough(biRegModelRawOut$output$status$code, 0)

omxCheckCloseEnough(biRegModelRawOut$output$iterations, 30, 10)

biRegModelRawBoot <- mxBootstrap(biRegModelRawOut, 10)
omxCheckTrue(is.null(biRegModelRawBoot$output[["standardErrors"]]))
bq1 <- summary(biRegModelRawBoot)[["bootstrapQuantile"]]
omxCheckCloseEnough(cor(bq1[,2] - bq1[,1],
                        biRegModelRawOut$output$standardErrors), 1.0, .45)

biRegModelRawBoot <- mxBootstrap(biRegModelRawBoot)
bq2 <- summary(biRegModelRawBoot)[["bootstrapQuantile"]]
omxCheckCloseEnough(cor(bq2[,2] - bq2[,1],
                        biRegModelRawOut$output$standardErrors), 1.0, .01)
repl3 <- biRegModelRawBoot$compute$output$raw[3,]

biRegModelRawBoot <- mxBootstrap(biRegModelRawBoot, 10)
bq3 <- summary(biRegModelRawBoot)[["bootstrapQuantile"]]
omxCheckEquals(bq3, bq1)

biRegModelRawBoot3 <- mxBootstrap(biRegModelRawBoot, only=3)
omxCheckEquals(repl3, biRegModelRawBoot3$compute$output$raw)

# investigate replication 3
biRegModelRaw3 <- mxModel(
  biRegModelRaw,
  mxData(observed=cbind(multiData1,
                        weight=biRegModelRawBoot3$compute$output$weight[[1]]),
         type="raw", weight = "weight"))

biRegModelRaw3 <- mxRun(biRegModelRaw3)

omxCheckCloseEnough(coef(biRegModelRaw3),
                    c(repl3[,names(coef(biRegModelRaw3))]), 1e-5)
