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


# -----------------------------------------------------------------------------
# Program: SimpleRegression_PathRaw.R  
# Author: Ryne Estabrook
# Date: 2009.08.01 
#
# ModelType: Regression
# DataType: Continuous
# Field: None
#
# Purpose: 
#      Simple Regression model to estimate effect of independent on dependent variables
#      Path style model input - Raw data input
#
# RevisionHistory:
#      Hermine Maes -- 2009.10.08 updated & reformatted
#      Ross Gore -- 2011.06.06 added Model, Data & Field metadata
#      Hermine Maes -- 2014.11.02 piecewise specification
# -----------------------------------------------------------------------------

require(OpenMx)
# Load Library
# -----------------------------------------------------------------------------

data(myRegDataRaw)

names(myRegDataRaw)

myRegDataRaw <- myRegDataRaw[,c("x","y")]
# Prepare Data
# -----------------------------------------------------------------------------

# dataset
dataRaw      <- mxData( observed=myRegDataRaw, type="raw" )
# variance paths
varPaths     <- mxPath( from=c("x", "y"), arrows=2, 
                        free=TRUE, values = c(1, 1), labels=c("varx", "residual") )
# regression weights
regPaths     <- mxPath( from="x", to="y", arrows=1, 
                        free=TRUE, values=1, labels="beta1" ) 
# means and intercepts
means        <- mxPath( from="one", to=c("x", "y"), arrows=1, 
                        free=TRUE, values=c(1, 1), labels=c("meanx", "beta0") )
    
uniRegModel  <- mxModel(model="Simple Regression Path Specification", type="RAM", 
                        dataRaw, manifestVars=c("x", "y"), varPaths, regPaths, means)                        
# Create an MxModel object
# -----------------------------------------------------------------------------
      
uniRegFit <- mxRun(uniRegModel, suppressWarnings=TRUE)

summary(uniRegFit)
uniRegFit$output


omxCheckCloseEnough(uniRegFit$output$estimate[["beta0"]], 2.5478, 0.001)
omxCheckCloseEnough(uniRegFit$output$estimate[["beta1"]], 0.4831, 0.001)
omxCheckCloseEnough(uniRegFit$output$estimate[["residual"]], 0.6652, 0.001)
omxCheckCloseEnough(uniRegFit$output$estimate[["meanx"]], 0.0542, 0.001)
omxCheckCloseEnough(uniRegFit$output$estimate[["varx"]], 1.1053, 0.001)
# Compare OpenMx results to Mx results 
# -----------------------------------------------------------------------------
