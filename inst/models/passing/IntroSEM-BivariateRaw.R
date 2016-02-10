# ---------------------------------------------------------------------
# Program: BivariateRaw-OpenMx100214.R
#  Author: Steven M. Boker
#    Date: Sun Feb 14 13:27:16 EST 2010
#
# This program fits a FIML bivariate model to the 
#     multiData simulated data.
#
#
# ---------------------------------------------------------------------
# Revision History
#    -- Sun Feb 14 13:27:13 EST 2010
#      Created BivariateRaw-OpenMx100214.R.
#    -- Sun Aug 29 2010
#      Formatted for OpenMx Test Suite
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
    mxData(observed=multiData1, type="raw")
    )

omxCheckError(mxRun(biRegModelRaw, suppressWarnings=TRUE),
"The job for model 'FIML Bivariate Regression of y on x1 and x2' exited abnormally with the error message: FIML Bivariate Regression of y on x1 and x2.expectation: raw data observed but no expected means vector was provided. Add something like mxPath(from = 'one', to = manifests) to your model.")

biRegModelRaw <- mxModel(
  biRegModelRaw,
  mxPath(from="one", to=manifests, 
         arrows=1, free=TRUE, values=.1, 
         labels=c("MeanX1", "MeanX2", "MeanY")))

biRegModelRawOut <- mxRun(biRegModelRaw)

summary(biRegModelRawOut)


# ----------------------------------
# check for correct values

expectVal <- c(0.4479, 0.4328, 1.1364, 0.5811, 1.5556, 1.412, 0.9849, 
1.9741, 2.5529)

expectSE <- c(0.0554, 0.0474, 0.0719, 0.0649, 0.0984, 0.0893, 0.0477, 0.0558, 
0.1004)

expectMin <- 4608.207

omxCheckCloseEnough(expectVal, biRegModelRawOut$output$estimate, 0.001)

omxCheckCloseEnough(expectSE, 
    as.vector(biRegModelRawOut$output[['standardErrors']]), 0.001)

omxCheckCloseEnough(expectMin, biRegModelRawOut$output$minimum, 0.001)

omxCheckCloseEnough(biRegModelRawOut$output$status$code, 0)

omxCheckCloseEnough(biRegModelRawOut$output$iterations, 30, 10)
