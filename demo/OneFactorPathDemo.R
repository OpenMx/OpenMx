# ---------------------------------------------------------------------
# Program: OneFactorPathDemo.R
#  Author: Steve Boker
#    Date: Thu Jul 30 13:33:08 EDT 2009
#
# This program is the OpenMx one factor path model demo for the front page
#
#
# ---------------------------------------------------------------------
# Revision History
#    -- Thu Jul 30 13:33:11 EDT 2009
#      Created OneFactorPathDemo.R.
#
# ---------------------------------------------------------------------

require(OpenMx)
demoData <- data("demoOneFactor.csv", header=T)
manifests <- names(demoData)
latents <- c("G")
factorModel <- mxModel("One Factor", type="RAM",
    manifestVars = manifests, latentVars = latents,
    mxPath(from=latents, to=manifests),
    mxPath(from=manifests, arrows=2),
    mxPath(from=latents, arrows=2, free=F, values=1.0),
    mxData(cov(demoData), type="cov", numObs=500))
summary(mxRun(factorModel))

