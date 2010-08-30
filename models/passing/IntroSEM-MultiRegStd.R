# ---------------------------------------------------------------------
# Program: MultiRegStd-OpenMx100214.R
#  Author: Steven M. Boker
#    Date: Sun Feb 14 12:25:16 EST 2010
#
# This program fits a multiple regression model to the 
#    multiData simulated data.
#
#
# ---------------------------------------------------------------------
# Revision History
#    -- Sun Feb 14 12:25:21 EST 2010
#      Created MultiRegStd-OpenMx100214.R.
#
# ---------------------------------------------------------------------

# ----------------------------------
# Read libraries and set options.

require(OpenMx)

# ----------------------------------
# Read the data and print descriptive statistics.

data(multiData1)

# ----------------------------------
# Build an OpenMx multiple regression model using y and x1

predictors <- c("x1", "x2", "x3", "x4")
outcomes <- c("y")
manifests <- names(multiData1)
multiData1Cov <- cov(multiData1)

multiRegModel <- mxModel("Multiple Regression of y on x1, x2, x3, and x4",
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
           arrows=2, all=TRUE,
           free=TRUE, values=.2),
    mxPath(from=manifests, 
           arrows=2, 
           free=TRUE, values=.8, 
           labels=c("VarX1", "VarX2", "VarX3", "VarX4", "VarE")),
    mxData(observed=multiData1Cov, type="cov", numObs=500)
    )

multiRegModelOut <- mxRun(multiRegModel)

summary(multiRegModelOut)


# ----------------------------------
# check for correct values

expectVal <- c(0.04427, 0.306983, 0.398635, 0.471855, 1.138705, 
0.582269, 1.558721, 0.636251, 0.566045, 2.106981, 0.591722, 0.435226, 
0.651973, 2.558098, 0.534445)

expectSE <- c(0.037002, 0.029577, 0.025352, 0.022138, 0.072093, 0.065091, 
0.098683, 0.074968, 0.084998, 0.133396, 0.080868, 0.091494, 0.107955, 
0.161951, 0.033835)

expectMin <- 3019.388

omxCheckCloseEnough(expectVal, multiRegModelOut@output$estimate, 0.001)

omxCheckCloseEnough(expectSE, 
    as.vector(multiRegModelOut@output$standardError), 0.001)

omxCheckCloseEnough(expectMin, multiRegModelOut@output$minimum, 0.001)


    

