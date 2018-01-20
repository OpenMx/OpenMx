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
# Program: RAM-3Factor-12Indicators.R
#  Author: Steven M. Boker
#    Date: Fri Jul 30 13:45:12 EDT 2010
#
# This program is a factor model using standard RAM.
#
# ---------------------------------------------------------------------
# Revision History
#    -- Fri Jul 30 13:45:12 EDT 2010
#      Created RAM-3Factor-12Indicators.R.
#
# ---------------------------------------------------------------------

# ----------------------------------
# Read libraries and set options.

library(OpenMx)

options(width=100)
set.seed(10)

# ---------------------------------------------------------------------
# Data for factor model.

numberSubjects <- 1000
numberFactors <- 3
numberIndPerFactor <- 8
numberIndicators <- numberIndPerFactor*numberFactors # must be a multiple of numberFactors

XMatrix <- matrix(rnorm(numberSubjects*numberFactors, mean=0, sd=1), numberSubjects, numberFactors)

tLoadings <- c(1, seq(.5, .9, length.out=(numberIndPerFactor-1)), rep(0, numberIndPerFactor*2),
  rep(0, numberIndPerFactor*1), 1, seq(.5, .9, length.out=(numberIndPerFactor-1)), rep(0, numberIndPerFactor*1),
  rep(0, numberIndPerFactor*2), 1, seq(.5, .9, length.out=(numberIndPerFactor-1)))
BMatrix <- matrix(tLoadings, numberFactors, numberIndicators, byrow=TRUE)
UMatrix <- matrix(rnorm(numberSubjects*numberIndicators, mean=0, sd=1), numberSubjects, numberIndicators)
YMatrix <- XMatrix %*% BMatrix + UMatrix

cor(XMatrix)

dimnames(YMatrix) <- list(NULL, paste("X", 1:numberIndicators, sep=""))

round(cor(YMatrix), 3)
round(cov(YMatrix), 3)

indicators <- paste("X", 1:numberIndicators, sep="")
totalVars <- numberIndicators + numberFactors

# ----------------------------------
# Build an orthogonal simple structure factor model

latents <- paste("F", 1:numberFactors, sep="")

uniqueLabels <- paste("U_", indicators, sep="")
meanLabels <- paste("M_", latents, sep="")
factorVarLabels <- paste("Var_", latents, sep="")

latents1 <- latents[1]
indicators1 <- indicators[1:numberIndPerFactor]
loadingLabels1 <- paste("b_F1", indicators[1:numberIndPerFactor], sep="") 
latents2 <- latents[2]
indicators2 <- indicators[numberIndPerFactor+(1:numberIndPerFactor)]
loadingLabels2 <- paste("b_F2", indicators[numberIndPerFactor+(1:numberIndPerFactor)], sep="") 
latents3 <- latents[3]
indicators3 <- indicators[(2*numberIndPerFactor)+(1:numberIndPerFactor)]
loadingLabels3 <- paste("b_F3", indicators[(2*numberIndPerFactor)+(1:numberIndPerFactor)], sep="") 

threeFactorOrthogonal <- mxModel("threeFactorOrthogonal",
    type="RAM",
    manifestVars=c(indicators),
    latentVars=c(latents,"dummy1"),
    mxPath(from=latents1, to=indicators1, 
           arrows=1, connect="all.pairs",
           free=TRUE, values=.2, ubound=5,
           labels=loadingLabels1),
    mxPath(from=latents2, to=indicators2, 
           arrows=1, connect="all.pairs",
           free=TRUE, values=.2, 
           labels=loadingLabels2),
    mxPath(from=latents3, to=indicators3, 
           arrows=1, connect="all.pairs",
           free=TRUE, values=.2, 
           labels=loadingLabels3),
    mxPath(from=latents1, to=indicators1[1], 
           arrows=1, 
           free=FALSE, values=1),
    mxPath(from=latents2, to=indicators2[1], 
           arrows=1, 
           free=FALSE, values=1),
    mxPath(from=latents3, to=indicators3[1], 
           arrows=1, 
           free=FALSE, values=1),
    mxPath(from=indicators, 
           arrows=2, 
           free=TRUE, values=.2, 
           labels=uniqueLabels),
    mxPath(from=latents,
           arrows=2, 
           free=TRUE, values=.8, lbound=1e-5,
           labels=factorVarLabels),
    mxPath(from="one", to=indicators, 
           arrows=1, free=FALSE, values=0),
    mxPath(from="one", to=c(latents), 
           arrows=1, free=TRUE, values=.1, 
           labels=meanLabels),
    mxCI(c('A', 'S'), boundAdj=FALSE),
    mxData(observed=cov(YMatrix), means=apply(YMatrix, 2, mean), 
	numObs=nrow(YMatrix), type="cov")
)

threeFactorOrthogonalOut <- mxRun(threeFactorOrthogonal)
omxCheckCloseEnough(threeFactorOrthogonalOut$output$fit, 29350.18, .1)

threeFactorCI <- omxParallelCI(threeFactorOrthogonalOut)

#cat(deparse(fivenum(threeFactorCI$output$confidenceIntervals[,'lbound'])))
omxCheckCloseEnough(fivenum(threeFactorCI$output$confidenceIntervals[,'lbound']),
                    c(0.389, 0.631, 0.837, 0.935, 0.987), .01)
omxCheckCloseEnough(fivenum(threeFactorCI$output$confidenceIntervals[,'ubound']),
                    c(0.545, 0.804, 1.04, 1.145, 1.344), .01)
