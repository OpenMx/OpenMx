# -----------------------------------------------------------------------
# Program: MultivariateRegression_PathRaw.R  
#  Author: Ryne Estabrook
#    Date: 08 01 2009 
#
# Multivariate Regression model to estimate effect of independent on dependent variables
# Path style model input - Raw data input
#
# Revision History
#   Hermine Maes -- 10 08 2009 updated & reformatted
# -----------------------------------------------------------------------

require(OpenMx)

#Prepare Data
# -----------------------------------------------------------------------
data(myRegDataRaw)

#Create an MxModel object
# -----------------------------------------------------------------------
multivariateRegModel <- mxModel("MultiVariate Regression -- Path Specification", 
    type="RAM",
    mxData(
        observed=myRegDataRaw, 
        type="raw"
    ),
    manifestVars=c("w", "x", "y", "z"),
    # variance paths
    mxPath(
        from=c("w", "x", "y", "z"), 
        arrows=2,
        free=TRUE, 
        values = c(1, 1, 1),
        labels=c("residualw", "varx", "residualy", "varz")
    ),
    # covariance of x and z
    mxPath(
        from="x",
        to="z",
        arrows=2,
        free=TRUE,
        values=0.5,
        labels="covxz"
    ), 
    # regression weights for y
    mxPath(
        from=c("x","z"),
        to="y",
        arrows=1,
        free=TRUE,
        values=1,
        labels=c("betayx","betayz")
    ), 
    # regression weights for w
    mxPath(
        from=c("x","z"),
        to="w",
        arrows=1,
        free=TRUE,
        values=1,
        labels=c("betawx","betawz")
    ), 
    # means and intercepts
    mxPath(
        from="one",
        to=c("w", "x", "y", "z"),
        arrows=1,
        free=TRUE,
        values=c(1, 1),
        labels=c("betaw", "meanx", "betay", "meanz")
    )
) # close model
  
multivariateRegFit <- mxRun(multivariateRegModel)

summary(multivariateRegFit)  
multivariateRegFit@output

#Compare OpenMx results to Mx results 
# -----------------------------------------------------------------------
omxCheckCloseEnough(multivariateRegFit@output$estimate[["betay"]], 1.6332, 0.001)
omxCheckCloseEnough(multivariateRegFit@output$estimate[["betayx"]], 0.4246, 0.001)
omxCheckCloseEnough(multivariateRegFit@output$estimate[["betayz"]], 0.2260, 0.001)
omxCheckCloseEnough(multivariateRegFit@output$estimate[["residualy"]], 0.6267, 0.001)
omxCheckCloseEnough(multivariateRegFit@output$estimate[["betaw"]], 0.5139, 0.001)
omxCheckCloseEnough(multivariateRegFit@output$estimate[["betawx"]], -0.2310, 0.001)
omxCheckCloseEnough(multivariateRegFit@output$estimate[["betawz"]], 0.5122, 0.001)
omxCheckCloseEnough(multivariateRegFit@output$estimate[["residualw"]], 0.5914, 0.001)
omxCheckCloseEnough(multivariateRegFit@output$estimate[["varx"]], 1.1053, 0.001)
omxCheckCloseEnough(multivariateRegFit@output$estimate[["varz"]], 0.8275, 0.001)
omxCheckCloseEnough(multivariateRegFit@output$estimate[["covxz"]], 0.2862, 0.001)
omxCheckCloseEnough(multivariateRegFit@output$estimate[["meanx"]], 0.0542, 0.001)
omxCheckCloseEnough(multivariateRegFit@output$estimate[["meanz"]], 4.0611, 0.001)
