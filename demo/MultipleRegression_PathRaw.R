#
#   Copyright 2007-2018 The OpenMx Project
#
#   Licensed under the Apache License, Version 2.0 (the "License");
#   you may not use this file except in compliance with the License.
#   You may obtain a copy of the License at
# 
#        http://www.apache.org/licenses/LICENSE-2.0
# 
#   Unless required by applicable law or agreed to in writing, software
#   distributed under the License is distributed on an "AS IS" BASIS,
#   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#   See the License for the specific language governing permissions and
#   limitations under the License.

require(OpenMx)

# -----------------------------------------------------------------------------
# Program: MultipleRegression_PathRaw.R  
# Author: Ryne Estabrook
# Date: 2009.08.01 
#
# ModelType: Regression
# DataType: Continuous
# Field: None
#
# Purpose:
#      Multiple Regression model to estimate effect of independent on dependent variables
#      Path style model input - Raw data input
#
# RevisionHistory:
#      Hermine Maes -- 2009.10.08 updated & reformatted
#      Ross Gore -- 2011.06.15 added Model, Data & Field metadata
#      Hermine Maes -- 2014.11.02 piecewise specification
# -----------------------------------------------------------------------------

require(OpenMx)
# Load Library
# -----------------------------------------------------------------------------

data(myRegDataRaw)

myRegDataRaw<-myRegDataRaw[,c("x","y","z")]
# Prepare Data
# -----------------------------------------------------------------------------

# dataset
dataRaw      <- mxData( observed=myRegDataRaw, type="raw" )
# variance paths      
varPaths     <- mxPath( from=c("x","y","z"),  arrows=2, 
                        free=TRUE, values = c(1,1,1), labels=c("varx","res","varz") )
# covariance of x and z
covPaths     <- mxPath( from="x", to="z", arrows=2, 
                        free=TRUE, values=0.5, labels="covxz" )
# regression weights
regPaths     <- mxPath( from=c("x","z"), to="y", arrows=1, 
                        free=TRUE, values=1, labels=c("betax","betaz") )
# means and intercepts
means        <- mxPath( from="one", to=c("x","y","z"), arrows=1, 
                        free=TRUE, values=c(1,1), labels=c("meanx","beta0","meanz") )

multiRegModel <- mxModel("Multiple Regression Path Specification", type="RAM",
                        dataRaw, manifestVars=c("x","y","z"), 
                        varPaths, covPaths, regPaths, means)
# Create an MxModel object
# -----------------------------------------------------------------------------
      
multiRegFit<-mxRun(multiRegModel)

summary(multiRegFit)
multiRegFit$output


omxCheckCloseEnough(multiRegFit$output$estimate[["beta0"]], 1.6332, 0.001)
omxCheckCloseEnough(multiRegFit$output$estimate[["betax"]], 0.4246, 0.001)
omxCheckCloseEnough(multiRegFit$output$estimate[["betaz"]], 0.2260, 0.001)
omxCheckCloseEnough(multiRegFit$output$estimate[["res"]], 0.6267, 0.001)
omxCheckCloseEnough(multiRegFit$output$estimate[["varx"]], 1.1053, 0.001)
omxCheckCloseEnough(multiRegFit$output$estimate[["varz"]], 0.8275, 0.001)
omxCheckCloseEnough(multiRegFit$output$estimate[["covxz"]], 0.2862, 0.001)
omxCheckCloseEnough(multiRegFit$output$estimate[["meanx"]], 0.0542, 0.001)
omxCheckCloseEnough(multiRegFit$output$estimate[["meanz"]], 4.0611, 0.001)
# Compare OpenMx results to Mx results 
# -----------------------------------------------------------------------------