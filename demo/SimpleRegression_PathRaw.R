# -----------------------------------------------------------------------
# Program: SimpleRegression_PathRaw.R  
#  Author: Ryne Estabrook
#    Date: 08 01 2009 
#
# Simple Regression model to estimate effect of independent on dependent variables
# Path style model input - Raw data input
#
# Revision History
#   Hermine Maes -- 10 08 2009 updated & reformatted
# -----------------------------------------------------------------------

require(OpenMx)

#Prepare Data
# -----------------------------------------------------------------------
data(myRegDataRaw)

names(myRegDataRaw)

myRegDataRaw <- myRegDataRaw[,c("x","y")]

#Create an MxModel object
# -----------------------------------------------------------------------
uniRegModel <- mxModel("Simple Regression -- Path Specification", 
      type="RAM",
      mxData(
      		observed=myRegDataRaw, 
      		type="raw"
      ),
      manifestVars=c("x","y"),
      # variances
      mxPath(from=c("x","y"), 
            arrows=2,
            free=TRUE, 
            values = c(1,1),
            labels=c("varx","residual")
      ), 
      # regression weight
      mxPath(from="x",
            to="y",
            arrows=1,
            free=TRUE,
            values=1,
            label="beta1"
      ), 
      # means
      mxPath(from="one",
            to=c("x","y"),
            arrows=1,
            free=TRUE,
            values=c(1,1),
            labels=c("meanx","beta0")
      ) 
) # close model
      
uniRegFit <- mxRun(uniRegModel)

summary(uniRegFit)
uniRegFit@output

#Compare OpenMx results to Mx results 
# -----------------------------------------------------------------------
omxCheckCloseEnough(uniRegFit@output$estimate[["beta0"]], 2.5478, 0.001)
omxCheckCloseEnough(uniRegFit@output$estimate[["beta1"]], 0.4831, 0.001)
omxCheckCloseEnough(uniRegFit@output$estimate[["residual"]], 0.6652, 0.001)
omxCheckCloseEnough(uniRegFit@output$estimate[["meanx"]], 0.0542, 0.001)
omxCheckCloseEnough(uniRegFit@output$estimate[["varx"]], 1.1053, 0.001)

