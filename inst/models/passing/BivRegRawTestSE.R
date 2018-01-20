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


# ----------------------------------
# Read libraries and set options.

options(width=80)

require(OpenMx)

# ----------------------------------
# Read the data and print descriptive statistics.

multiData1 <- suppressWarnings(try(read.csv("models/passing/data/multiData.csv"), silent=TRUE))
if (is(multiData1, "try-error")) multiData1 <- read.csv("data/multiData.csv")


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
    mxPath(from="one", to=manifests, 
           arrows=1, free=TRUE, values=.1, 
           labels=c("MeanX1", "MeanX2", "MeanY")),
    mxData(observed=multiData1, type="raw")
    )
biRegModelRaw <- mxOption(biRegModelRaw, "Standard Errors", "Yes")

biRegModelRawOut <- mxRun(biRegModelRaw)

summary(biRegModelRawOut)

lmOut <- lm(y~x1+x2, data=multiData1)
lmOutSummary <- summary(lmOut)

omxCheckCloseEnough(biRegModelRawOut$output$estimate[1:2], 
    lmOutSummary$coef[2:3,1],
    0.001)

omxCheckCloseEnough(biRegModelRawOut$output$standardErrors[1:2], 
    lmOutSummary$coef[2:3,2],
    0.001)

