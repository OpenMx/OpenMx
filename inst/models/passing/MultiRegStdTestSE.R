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

multiData1 <- read.csv("data/multiData.csv")


# ----------------------------------
# Build an OpenMx multiple regression model using y and x1

predictors <- c("x1", "x2", "x3", "x4")
outcomes <- c("y")
manifests <- names(multiData1)
multiData1Cov <- cov(multiData1)

multiRegModel <- mxModel("Multiple Regression of y on x1, x2, x3, x4",
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
multiRegModel <- mxOption(multiRegModel, "Standard Errors", "Yes")

multiRegModelOut <- mxRun(multiRegModel)

summary(multiRegModelOut)

lmOut <- lm(y~x1+x2+x3+x4, data=multiData1)
lmOutSummary <- summary(lmOut)

omxCheckCloseEnough(multiRegModelOut$output$estimate[1:4], 
    lmOutSummary$coef[2:5,1],
    0.001)

omxCheckCloseEnough(multiRegModelOut$output$standardErrors[1:4], 
    lmOutSummary$coef[2:5,2],
    0.001)
