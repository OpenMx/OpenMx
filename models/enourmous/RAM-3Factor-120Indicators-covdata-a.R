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
numberIndPerFactor <- 40
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

YFrame <- data.frame(YMatrix)

round(cor(YFrame), 3)
round(cov(YFrame), 3)

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
           arrows=1, all=TRUE, 
           free=TRUE, values=.2, 
           labels=loadingLabels1),
    mxPath(from=latents2, to=indicators2, 
           arrows=1, all=TRUE, 
           free=TRUE, values=.2, 
           labels=loadingLabels2),
    mxPath(from=latents3, to=indicators3, 
           arrows=1, all=TRUE, 
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
           free=TRUE, values=.8, 
           labels=factorVarLabels),
    mxPath(from="one", to=indicators, 
           arrows=1, free=FALSE, values=0),
    mxPath(from="one", to=c(latents), 
           arrows=1, free=TRUE, values=.1, 
           labels=meanLabels),
    mxData(observed=cov(YFrame), means=mean(YFrame), 
	numObs=nrow(YFrame), type="cov")
    )

threeFactorOrthogonalOut <- mxRun(threeFactorOrthogonal)
summary(threeFactorOrthogonalOut)

