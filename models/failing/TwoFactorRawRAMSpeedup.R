# ---------------------------------------------------------------------
# Program: twoFactorRawRAMSpeedup.R
#  Author: Steven M. Boker
#    Date: Thu Mar 4 20:23:00 EST 2010
#
# This program fits a FIML single factor model to simulated data.
#    using (a) standard RAM, (b) Geometric Series RAM.
#
# ---------------------------------------------------------------------
# Revision History
#    -- Thu Mar 4 20:23:04 EST 2010
#      Created twoFactorRawRAMSpeedup.R.
#
# ---------------------------------------------------------------------

# ----------------------------------
# Read libraries and set options.

options(width=100)

require(psych)
require(OpenMx)

# ----------------------------------
# Read the data and print descriptive statistics.

numberSubjects <- 200
numberIndicators <- 20
numberFactors <- 2

BMatrix <- matrix(seq(-1,1, length.out=numberIndicators), numberFactors, numberIndicators, byrow=TRUE)
XMatrix <- matrix(rnorm(numberSubjects*numberFactors, mean=0, sd=1), numberSubjects, numberFactors)
UMatrix <- matrix(rnorm(numberSubjects*numberIndicators, mean=0, sd=1), numberSubjects, numberIndicators)
YMatrix <- XMatrix %*% BMatrix + UMatrix

indicators <- paste("Y", 1:numberIndicators, sep="")
dimnames(YMatrix) <- list(NULL, indicators)
YFrame <- data.frame(YMatrix)

describe(YFrame)

# ----------------------------------
# Build an Old-style RAM OpenMx single factor FIML model with fixed variance

latents <- paste("F", 1:numberFactors, sep="")
loadingLabels <- paste("b_", indicators, sep="")
uniqueLabels <- paste("U_", indicators, sep="")
meanLabels <- paste("M_", indicators, sep="")
factorVarLabels <- paste("Var_", latents, sep="")

twoFactorRawOld <- mxModel("twoFactorOldRAM",
    type="RAM",
    manifestVars=indicators,
    latentVars=latents,
    mxPath(from=latents, to=indicators, 
           arrows=1, all=TRUE, 
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
    mxData(observed=YFrame, type="raw")
    )

twoFactorRawOldOut <- mxRun(twoFactorRawOld)

summary(twoFactorRawOldOut)

# ----------------------------------
# Build an New-style RAM OpenMx single factor FIML model with fixed variance


twoFactorRawNewRAM <- mxModel("twoFactorNewRAM",
    type="RAM",
    manifestVars=indicators,
    latentVars=latents,
    mxPath(from=latents, to=indicators, 
           arrows=1, all=TRUE, 
           free=TRUE, values=.2, 
           labels=loadingLabels),
    mxPath(from=indicators, 
           arrows=2, 
           free=TRUE, values=.8, 
           labels=uniqueLabels),
    mxPath(from=latents,
           arrows=2, 
           free=FALSE, values=1, 
           labels=factorVarLabels))
twoFactorRawNew <- mxModel(twoFactorRawNewRAM,
    type="raw",
    mxMatrix("Iden", length(indicators)+length(latents), name="I"),
    mxAlgebra(F %*% (I + A), name="E"),
    mxAlgebra(E %*% S %*% t(E), name="R", 
        dimnames = list(indicators, indicators)),
    mxMatrix("Full", nrow=1, ncol=length(indicators),
        values=0.1,
        free=TRUE,
        labels=meanLabels,
        dimnames=list(NULL, indicators),
        name="M"
    ),
    mxFIMLObjective(covariance="R", means="M"),
    mxData(observed=YFrame, type="raw")
    )

twoFactorRawNewOut <- mxRun(twoFactorRawNew)

summary(twoFactorRawNewOut)

summary(twoFactorRawOldOut)$parameters[,5:6] - summary(twoFactorRawNewOut)$parameters[,5:6] 

summary(twoFactorRawOldOut)$frontendTime - summary(twoFactorRawNewOut)$frontendTime

summary(twoFactorRawOldOut)$backendTime - summary(twoFactorRawNewOut)$backendTime


