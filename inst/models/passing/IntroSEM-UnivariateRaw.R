# ---------------------------------------------------------------------
# Program: UnivariateRaw-OpenMx100214.R
#  Author: Steven M. Boker
#    Date: Sun Feb 14 13:18:47 EST 2010
#
# This program fits a FIML univariate model to the 
#     multiData simulated data.
#
#
# ---------------------------------------------------------------------
# Revision History
#    -- Sun Feb 14 13:18:51 EST 2010
#      Created UnivariateRaw-OpenMx100214.R.
#
# ---------------------------------------------------------------------

# ----------------------------------
# Read libraries and set options.

require(OpenMx)

# ----------------------------------
# Read the data and print descriptive statistics.

data(multiData1)

# ----------------------------------
# Build an OpenMx univariate regression model using y and x1

manifests <- c("x1", "y")

uniRegModelRaw <- mxModel("FIML Univariate Regression of y on x1",
    type="RAM",
    manifestVars=manifests,
    mxPath(from="x1", to="y", arrows=1, 
           free=TRUE, values=.2, labels="b1"),
    mxPath(from=manifests, 
           arrows=2, free=TRUE, values=.8, 
           labels=c("VarX1", "VarE")),
    mxPath(from="one", to=manifests, 
           arrows=1, free=TRUE, values=.1, 
           labels=c("MeanX1", "MeanY")),
    mxData(observed=multiData1, type="raw")
    )

uniRegModelRawOut <- mxRun(uniRegModelRaw, suppressWarnings=TRUE)

summary(uniRegModelRawOut)


#---------------------
# check values: uniRegModelRawOut

expectVal <- c(0.669179, 1.13643, 1.647629, 0.984894, 3.189368)

expectSE <-c(0.053849, 0.071873, 0.104204, 0.047674, 0.078154)

expectMin <- 3151.492

omxCheckCloseEnough(expectVal, uniRegModelRawOut$output$estimate, 0.001)

omxCheckCloseEnough(expectSE, 
    as.vector(uniRegModelRawOut$output[['standardErrors']]), 0.001)

omxCheckCloseEnough(expectMin, uniRegModelRawOut$output$minimum, 0.001)




    

