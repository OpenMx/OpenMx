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
# Program: MultivariateRegression_PathRaw.R  
# Author: Ryne Estabrook
# Date: 2009.08.01 
#
# ModelType: Regression
# DataType: Continuous
# Field: None
#
# Purpose: 
#      Multivariate Regression model to estimate effect of independent on dependent variables
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
# Prepare Data
# -----------------------------------------------------------------------------

# dataset
dataRaw      <- mxData( observed=myRegDataRaw, type="raw" )
# variance paths
varPaths     <- mxPath( from=c("w","x","y","z"), arrows=2, 
                        free=TRUE, values=1, 
                        labels=c("residualw","varx","residualy","varz") )
# covariance of x and z
covPaths     <- mxPath( from="x", to="z", arrows=2, 
                        free=TRUE, values=0.5, labels="covxz" ) 
# regression weights for y
regPathsY    <- mxPath( from=c("x","z"), to="y", arrows=1, 
                        free=TRUE, values=1, labels=c("betayx","betayz") ) 
# regression weights for w
regPathsW    <- mxPath( from=c("x","z"), to="w", arrows=1, 
                        free=TRUE, values=1, labels=c("betawx","betawz") ) 
# means and intercepts
means        <- mxPath( from="one", to=c("w","x","y","z"), arrows=1, 
                        free=TRUE, values=c(1, 1), 
                        labels=c("betaw","meanx","betay","meanz") )

multivariateRegModel <- mxModel("MultiVariate Regression Path Specification", 
                        type="RAM", dataRaw, manifestVars=c("w","x","y","z"),
                        varPaths, covPaths, regPathsY, regPathsW, means )
# Create an MxModel object
# -----------------------------------------------------------------------------
  
multivariateRegFit <- mxRun(multivariateRegModel)

summary(multivariateRegFit)  
multivariateRegFit$output


omxCheckCloseEnough(multivariateRegFit$output$estimate[["betay"]], 1.6332, 0.001)
omxCheckCloseEnough(multivariateRegFit$output$estimate[["betayx"]], 0.4246, 0.001)
omxCheckCloseEnough(multivariateRegFit$output$estimate[["betayz"]], 0.2260, 0.001)
omxCheckCloseEnough(multivariateRegFit$output$estimate[["residualy"]], 0.6267, 0.001)
omxCheckCloseEnough(multivariateRegFit$output$estimate[["betaw"]], 0.5139, 0.001)
omxCheckCloseEnough(multivariateRegFit$output$estimate[["betawx"]], -0.2310, 0.001)
omxCheckCloseEnough(multivariateRegFit$output$estimate[["betawz"]], 0.5122, 0.001)
omxCheckCloseEnough(multivariateRegFit$output$estimate[["residualw"]], 0.5914, 0.001)
omxCheckCloseEnough(multivariateRegFit$output$estimate[["varx"]], 1.1053, 0.001)
omxCheckCloseEnough(multivariateRegFit$output$estimate[["varz"]], 0.8275, 0.001)
omxCheckCloseEnough(multivariateRegFit$output$estimate[["covxz"]], 0.2862, 0.001)
omxCheckCloseEnough(multivariateRegFit$output$estimate[["meanx"]], 0.0542, 0.001)
omxCheckCloseEnough(multivariateRegFit$output$estimate[["meanz"]], 4.0611, 0.001)
# Compare OpenMx results to Mx results 
# -----------------------------------------------------------------------------