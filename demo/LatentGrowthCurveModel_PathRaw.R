# -----------------------------------------------------------------------
# Program: LatentGrowthModel_PathRaw.R  
#  Author: Ryne Estabrook
#    Date: 08 01 2009 
#
# Latent Growth model to estimate means and (co)variances of slope and intercept
# Path style model input - Raw data input
#
# Revision History
#   Hermine Maes -- 10 08 2009 updated & reformatted
# -----------------------------------------------------------------------

require(OpenMx)

#Prepare Data
# -----------------------------------------------------------------------
data(myLongitudinalData)

#Create an MxModel object
# -----------------------------------------------------------------------
growthCurveModel <- mxModel("Linear Growth Curve Model -- Path Specification", 
    type="RAM",
    mxData(
    	observed=myLongitudinalData,
        type="raw"
    ),
    manifestVars=c("x1","x2","x3","x4","x5"),
    latentVars=c("intercept","slope"),
    # residual variances
    mxPath(
    	from=c("x1","x2","x3","x4","x5"), 
        arrows=2,
        free=TRUE, 
        values = c(1, 1, 1, 1, 1),
        labels=c("residual","residual","residual","residual","residual")
    ),
    # latent variances and covariance
    mxPath(
    	from=c("intercept","slope"), 
        arrows=2,
        all=TRUE,
        free=TRUE, 
        values=c(1, 1, 1, 1),
        labels=c("vari", "cov", "cov", "vars")
    ),
    # intercept loadings
    mxPath(
    	from="intercept",
        to=c("x1","x2","x3","x4","x5"),
        arrows=1,
        free=FALSE,
        values=c(1, 1, 1, 1, 1)
    ),
    # slope loadings
    mxPath(
    	from="slope",
        to=c("x1","x2","x3","x4","x5"),
        arrows=1,
        free=FALSE,
        values=c(0, 1, 2, 3, 4)
    ),
    # manifest means
    mxPath(from="one",
        to=c("x1", "x2", "x3", "x4", "x5"),
        arrows=1,
        free=FALSE,
        values=c(0, 0, 0, 0, 0)
    ),
    # latent means
    mxPath(from="one",
        to=c("intercept", "slope"),
        arrows=1,
        free=TRUE,
        values=c(1, 1),
        labels=c("meani", "means")
    )
) # close model
      
growthCurveFit <- mxRun(growthCurveModel)

summary(growthCurveFit)
growthCurveFit@output$estimate

#Compare OpenMx results to Mx results 
# -----------------------------------------------------------------------
omxCheckCloseEnough(growthCurveFit@output$estimate[["meani"]], 9.930, 0.01)
omxCheckCloseEnough(growthCurveFit@output$estimate[["means"]], 1.813, 0.01)
omxCheckCloseEnough(growthCurveFit@output$estimate[["vari"]], 3.878, 0.01)
omxCheckCloseEnough(growthCurveFit@output$estimate[["vars"]], 0.258, 0.01)
omxCheckCloseEnough(growthCurveFit@output$estimate[["cov"]], 0.460, 0.01)
omxCheckCloseEnough(growthCurveFit@output$estimate[["residual"]], 2.316, 0.01)
