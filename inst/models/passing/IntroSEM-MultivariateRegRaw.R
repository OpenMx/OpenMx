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


# ---------------------------------------------------------------------
# Program: MultivariateRegRaw-OpenMx100214.R
#  Author: Steven M. Boker
#    Date: Sun Feb 14 14:23:50 EST 2010
#
# This program fits a FIML multiple regression model to the 
#    multiData simulated data.
#
#
# ---------------------------------------------------------------------
# Revision History
#    -- Sun Feb 14 14:23:53 EST 2010
#      Created MultivariateRegRaw-OpenMx100214.R.
#
# ---------------------------------------------------------------------

# ----------------------------------
# Read libraries and set options.

library(OpenMx)

# ----------------------------------
# Read the data and print descriptive statistics.

data(multiData1)

# ----------------------------------
# Build an OpenMx multiple regression model using y and x1

predictors <- c("x2", "x3", "x4")
outcomes <- c("y", "x1")
manifests <- names(multiData1)

multivariateRegModel <- mxModel("FIML Multiple Regression of y on x1, x2, x3 and x4",
    type="RAM",
    manifestVars=manifests,
    mxPath(from=predictors, to=outcomes, 
#           arrows=1, all=TRUE,
           arrows=1, connect="all.pairs",
           free=TRUE, values=.2, 
           labels=c("b21", "b22", "b31", "b32", "b41", "b42")),
    mxPath(from="x1", to="y", 
           arrows=1,
           free=TRUE, values=.2, 
           labels="b11"),
    mxPath(from=predictors, to=predictors,
#           arrows=2, all=TRUE,
 			arrows=2, connect="unique.pairs",
           free=TRUE, values=.2),
    mxPath(from=manifests, to=manifests,
           arrows=2, 
           free=TRUE, values=.8, 
           labels=c("VarEx1", "VarX2", "VarX3", "VarX4", "VarEy")),
    mxPath(from="one", to=manifests, 
           arrows=1, free=TRUE, values=.1, 
           labels=c("MeanX1", "MeanX2", "MeanX3", "MeanX4", "MeanY")),
    mxData(observed=multiData1, type="raw")
    )

multivariateRegModelOut <- mxRun(multivariateRegModel, suppressWarnings=TRUE)

summary(multivariateRegModelOut)


# ----------------------------------
# check for correct values

expectVal <- c(0.04427, 0.26689, 0.306983, 0.187537, 0.398635, 0.138109, 
0.471855, 0.780697, 1.555604, 0.564914, 2.102769, 0.434357, 0.650671, 
2.552984, 0.533376, -0.678459, 1.974178, 3.073289, 4.055537, 
0.060049)

expectSE <- c(0.036963, 0.033694, 0.029546, 0.02947, 0.025327, 0.026034, 
0.022115, 0.049375, 0.098387, 0.084745, 0.13299, 0.091222, 0.107629, 
0.161465, 0.033733, 0.12533, 0.055778, 0.06485, 0.071456, 0.106583
)

expectMin <- 7615.126

omxCheckCloseEnough(expectMin, multivariateRegModelOut$output$minimum, 0.001)

omxCheckCloseEnough(expectVal, multivariateRegModelOut$output$estimate, 0.002)

omxCheckCloseEnough(expectSE, 
    as.vector(multivariateRegModelOut$output[['standardErrors']]), 0.001)

# ----------------------------

multiData1$x1o <- cut(multiData1$x1, breaks = 4, ordered_result = TRUE)

outcomes <- c('y', 'x1o')

r2 <- mxModel("regr", type="RAM",
              latentVars = predictors,
              manifestVars = outcomes,
              mxPath("one", predictors, free=FALSE, labels=paste0('data.', predictors)),
              mxPath(predictors, outcomes),
              mxPath("one", outcomes),
              mxPath(outcomes, arrows=2, values=1),
              mxThreshold(vars='x1o', nThresh=3),
              mxData(observed=multiData1, type="raw"),
              mxFitFunctionWLS(),
              mxComputeSequence(list(
                GD=mxComputeGradientDescent(),
                CK=mxComputeCheckpoint(toReturn = TRUE, vcov=TRUE, vcovWLS=TRUE, vcovFilter =
                                      c("regr.data.y.x2", "regr.data.y.x3",
                                        "regr.data.x1o.th1", "regr.data.x1o.th2")))))

r2 <- mxRun(r2)

l1 <- r2$compute$steps$CK$log

yv <- r2$data$observedStats$y.vcov
x1ov <- r2$data$observedStats$x1o.vcov

omxCheckCloseEnough(l1["regr.data.y.V(intercept):(intercept)"],
                    yv["(intercept)","(intercept)"])
omxCheckCloseEnough(l1["regr.data.y.Vx2:x2"],
                    yv["x2","x2"])
omxCheckCloseEnough(l1["regr.data.y.Vx3:x2"],
                    yv["x3","x2"])
omxCheckCloseEnough(l1["regr.data.y.Vx3:x3"],
                    yv["x3","x3"])
omxCheckCloseEnough(l1["regr.data.y.Vx4:x4"],
                    yv["x4","x4"])
omxCheckCloseEnough(l1["regr.data.x1o.Vth1:th1"],
                    x1ov["th1","th1"])
omxCheckCloseEnough(l1["regr.data.x1o.Vth2:th1"],
                    x1ov["th2","th1"])
omxCheckCloseEnough(l1["regr.data.x1o.Vth2:th2"],
                    x1ov["th2","th2"])
omxCheckCloseEnough(l1["regr.data.x1o.Vth3:th3"],
                    x1ov["th3","th3"])
omxCheckCloseEnough(l1["regr.data.x1o.Vx2:x2"],
                    x1ov["x2","x2"])
omxCheckCloseEnough(l1["regr.data.x1o.Vx3:x3"],
                    x1ov["x3","x3"])
omxCheckCloseEnough(l1["regr.data.x1o.Vx4:x4"],
                    x1ov["x4","x4"])
