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

numberSubjects <- 200
numberIndicators <- 10

BMatrix <- matrix(seq(-1,1, length.out=numberIndicators), 1, numberIndicators)
XMatrix <- matrix(rnorm(numberSubjects, mean=0, sd=1), numberSubjects, 1)
XMatrix <- (XMatrix-mean(XMatrix))/sd(XMatrix)
UMatrix <- matrix(rnorm(numberSubjects*numberIndicators, mean=0, sd=1), numberSubjects, numberIndicators)
YMatrix <- XMatrix %*% BMatrix + UMatrix

indicators <- paste("Y", 1:numberIndicators, sep="")
dimnames(YMatrix) <- list(NULL, indicators)
YFrame <- data.frame(YMatrix)

describe(YFrame)

# ----------------------------------
# Build an Old-style RAM OpenMx single factor FIML model with fixed variance

latents <- c("F1")
loadingLabels <- paste("b_", indicators, sep="")
uniqueLabels <- paste("U_", indicators, sep="")
meanLabels <- paste("M_", indicators, sep="")
factorVarLabels <- paste("Var_", latents, sep="")

oneFactorRawOld <- mxModel("OneFactorOldRAM",
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

oneFactorRawOldOut <- mxRun(oneFactorRawOld)

summary(oneFactorRawOldOut)

# ----------------------------------
# Build an New-style RAM OpenMx single factor FIML model with fixed variance

latents <- c("F1")
loadingLabels <- paste("b_", indicators, sep="")
uniqueLabels <- paste("U_", indicators, sep="")
meanLabels <- paste("M_", indicators, sep="")
factorVarLabels <- paste("Var_", latents, sep="")

oneFactorRawNewRAM <- mxModel("OneFactorNewRAM",
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
oneFactorRawNew <- mxModel(oneFactorRawNewRAM,
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

oneFactorRawNewOut <- mxRun(oneFactorRawNew)

summary(oneFactorRawNewOut)

summary(oneFactorRawOldOut)$parameters[,5:6] - summary(oneFactorRawNewOut)$parameters[,5:6] 

summary(oneFactorRawOldOut)$frontendTime - summary(oneFactorRawNewOut)$frontendTime

summary(oneFactorRawOldOut)$backendTime - summary(oneFactorRawNewOut)$backendTime


