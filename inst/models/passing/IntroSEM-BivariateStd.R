# ---------------------------------------------------------------------
# Program: BivariateStd-OpenMx100214.R
#  Author: Steven M. Boker
#    Date: Sun Feb 14 12:12:17 EST 2010
#
# This program fits a bivariate model to the multiData simulated data.
#
#
# ---------------------------------------------------------------------
# Revision History
#    -- Sun Feb 14 12:12:15 EST 2010
#      Created BivariateStd-OpenMx100214.R.
#
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
multiData1Cov <- cov(multiData1[,c(1,2,5)])

biRegModel <- mxModel("Bivariate Regression of y on x1 and x2",
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
    mxData(observed=multiData1Cov, type="cov", numObs=500)
    )

biRegModelOut <- mxRun(biRegModel, suppressWarnings=TRUE)

# ensure summary looks in model's runstate
biRegModelOut$compute$steps[["GD"]]$engine <- 'XYZ'

brmSum <- summary(biRegModelOut)
omxCheckCloseEnough(brmSum$CFI, 1, 1e-6)
omxCheckCloseEnough(brmSum$TLI, 1, 1e-6)
omxCheckCloseEnough(brmSum$RMSEA, 0, 1e-6)
omxCheckTrue(all(is.na(brmSum$RMSEACI)))

# ----------------------------------
# check for correct values

expectVal <- c(0.4479, 0.4327, 1.1387, 0.5823, 1.5587, 1.4148)

expectSE <- c(0.0555, 0.0474, 0.0721, 0.0651, 0.0987, 0.0896)

expectMin <- 1850.685

omxCheckCloseEnough(expectVal, biRegModelOut$output$estimate, 0.001)

omxCheckCloseEnough(expectSE, 
    as.vector(biRegModelOut$output[['standardErrors']]), 0.001)

omxCheckCloseEnough(expectMin, biRegModelOut$output$minimum, 0.001)

omxCheckEquals(brmSum$optimizerEngine, mxOption(NULL, "Default optimizer"))

