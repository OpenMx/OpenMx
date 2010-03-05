# ---------------------------------------------------------------------
# Program: OneFactorRawRAMSpeedup.R
#  Author: Steven M. Boker
#    Date: Thu Mar 4 17:53:39 EST 2010
#
# This program fits a FIML single factor model to simulated data.
#    using (a) standard RAM, (b) Geometric Series RAM.
#
# ---------------------------------------------------------------------
# Revision History
#    -- Thu Mar 4 17:53:43 EST 2010
#      Created OneFactorRawRAMSpeedup.R.
#
# ---------------------------------------------------------------------

# ----------------------------------
# Read libraries and set options.

options(width=100)

require(psych)
require(OpenMx)

# ----------------------------------
# Read the data and print descriptive statistics.

numberSubjects <- 1000
numberIndicators <- 40
numberFactors <- 3

BMatrix <- matrix(seq(-1,1, length.out=numberIndicators*numberFactors), numberFactors, numberIndicators, byrow=TRUE)
XMatrix <- matrix(rnorm(numberSubjects, mean=0, sd=1), numberSubjects, numberFactors)
UMatrix <- matrix(rnorm(numberSubjects*numberIndicators, mean=0, sd=1), numberSubjects, numberIndicators)
YMatrix <- XMatrix %*% BMatrix + UMatrix

indicators <- paste("Y", 1:numberIndicators, sep="")
dimnames(YMatrix) <- list(NULL, indicators)
YFrame <- data.frame(YMatrix)

describe(YFrame)

# ----------------------------------
# Build an Old-style RAM OpenMx single factor FIML model with fixed variance

latents <- paste("F", 1:numberFactors, sep="")
loadingLabels <- paste("b_F", rep(1:numberFactors, each=numberIndicators), rep(indicators, numberFactors), sep="") 
loadingLabels

uniqueLabels <- paste("U_", indicators, sep="")
meanLabels <- paste("M_", indicators, sep="")
factorVarLabels <- paste("Var_", latents, sep="")

totalVars <- numberIndicators + numberFactors

AMatrixFree <- matrix(FALSE, totalVars, totalVars)
AMatrixFree[1:numberIndicators,(numberIndicators+1):totalVars] <- TRUE
AMatrixVals <- matrix(0, totalVars, totalVars)
AMatrixVals[1:numberIndicators,(numberIndicators+1):totalVars] <- .2
AMatrixLabels <- matrix(NA, totalVars, totalVars)
AMatrixLabels[1:numberIndicators,(numberIndicators+1):totalVars] <- loadingLabels

SMatrixFree <- matrix(FALSE, totalVars, totalVars)
SMatrixFree[cbind(1:numberIndicators,1:numberIndicators)] <- TRUE
SMatrixVals <- matrix(0, totalVars, totalVars)
SMatrixVals[cbind(1:numberIndicators,1:numberIndicators)] <- .6
SMatrixVals[cbind((numberIndicators+1):totalVars,(numberIndicators+1):totalVars)] <- 1
SMatrixLabels <- matrix(NA, totalVars, totalVars)
SMatrixLabels[cbind(1:numberIndicators,1:numberIndicators)] <- uniqueLabels

oneFactorRawOld <- mxModel("oneFactorOldRAM",
    mxMatrix("Full", totalVars, totalVars, 
        labels=AMatrixLabels, 
        free=AMatrixFree, 
        values=AMatrixVals,
        name="A", 
        byrow=TRUE
    ),
    mxMatrix("Full", totalVars, totalVars, 
        labels=SMatrixLabels, 
        free=SMatrixFree, 
        values=SMatrixVals,
        name="S", 
        byrow=TRUE
    ),
    mxMatrix("Full", numberIndicators, totalVars, 
        free=FALSE, 
        values=cbind(diag(1,numberIndicators),matrix(0,numberIndicators,numberFactors)),
        name="F", 
        byrow=TRUE
    ),
    mxMatrix("Iden", totalVars, name="I"),
    mxAlgebra(F %*% solve(I-A) %*% S %*% t(solve(I-A)) %*% t(F), 
        name="R", 
        dimnames = list(indicators, indicators)
    ),
#    mxMatrix("Full", nrow=1, ncol=numberIndicators,
#        values=0.1,
#        free=TRUE,
#        labels=meanLabels,
#        dimnames=list(NULL, indicators),
#        name="M"
#    ),
#    mxFIMLObjective(covariance="R", means="M"),
#    mxData(YFrame, type="raw")
    mxMLObjective("R"),
    mxData(observed=cov(YFrame), type="cov", numObs=1000)
)

oneFactorRawOldOut <- mxRun(oneFactorRawOld)

tSummaryOld <- summary(oneFactorRawOldOut)
tSummaryOld

# ----------------------------------
# Build an New-style RAM OpenMx single factor FIML model with fixed variance


oneFactorRawNew <- mxModel(oneFactorRawOld,
    mxAlgebra(F %*% (I + A), name="E"),
    mxAlgebra(E %*% S %*% t(E), name="R", 
        dimnames = list(indicators, indicators))
    )

oneFactorRawNewOut <- mxRun(oneFactorRawNew)

tSummaryNew <- summary(oneFactorRawNewOut)
tSummaryNew

tSummaryOld$parameters[,5:6] - tSummaryNew$parameters[,5:6] 

tSummaryOld$frontendTime - tSummaryNew$frontendTime

tSummaryOld$backendTime - tSummaryNew$backendTime

as.double(tSummaryOld$backendTime - tSummaryNew$backendTime) / as.double(tSummaryOld$backendTime)

