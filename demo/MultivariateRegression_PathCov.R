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

# ------------------------------------------------------------------------------
# Program: MultivariateRegression_PathCov.R  
# Author: Ryne Estabrook
# Date: 2009.08.01 
#
# ModelType: Regression
# DataType: Continuous
# Field: None
#
# Purpose: 
#      Multivariate Regression model to estimate effect of 
#      independent on dependent variables
#      Path style model input - Covariance matrix data input
#
# RevisionHistory:
#      Hermine Maes -- 2009.10.08 updated & reformatted
#      Ross Gore -- 2011.06.15 added Model, Data & Field metadata
#      Hermine Maes -- 2014.11.02 piecewise specification
# -----------------------------------------------------------------------------

require(OpenMx)
# Load Library
# -----------------------------------------------------------------------------

myRegDataCov<-matrix(
	c(0.808,-0.110, 0.089, 0.361,
	 -0.110, 1.116, 0.539, 0.289,
	  0.089, 0.539, 0.933, 0.312,
	  0.361, 0.289, 0.312, 0.836),
	nrow=4,
	dimnames=list(
		c("w","x","y","z"),
		c("w","x","y","z"))
)
	
myRegDataMeans <- c(2.582, 0.054, 2.574, 4.061)
names(myRegDataMeans) <- c("w","x","y","z")
# Prepare Data
# -----------------------------------------------------------------------------

# dataset
dataCov      <- mxData( observed=myRegDataCov, type="cov", numObs=100, means=myRegDataMeans )
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
                        type="RAM", dataCov, manifestVars=c("w","x","y","z"),
                        varPaths, covPaths, regPathsY, regPathsW, means )
# Create an MxModel object
# -----------------------------------------------------------------------------
  
multivariateRegFit <- mxRun(multivariateRegModel)

summary(multivariateRegFit)
multivariateRegFit$output


omxCheckCloseEnough(multivariateRegFit$output$estimate[["betay"]], 1.6312, 0.001)
omxCheckCloseEnough(multivariateRegFit$output$estimate[["betayx"]], 0.4243, 0.001)
omxCheckCloseEnough(multivariateRegFit$output$estimate[["betayz"]], 0.2265, 0.001)
omxCheckCloseEnough(multivariateRegFit$output$estimate[["residualy"]], 0.6272, 0.001)
omxCheckCloseEnough(multivariateRegFit$output$estimate[["betaw"]], 0.5165, 0.001)
omxCheckCloseEnough(multivariateRegFit$output$estimate[["betawx"]], -0.2311, 0.001)
omxCheckCloseEnough(multivariateRegFit$output$estimate[["betawz"]], 0.5117, 0.001)
omxCheckCloseEnough(multivariateRegFit$output$estimate[["residualw"]], 0.5918, 0.001)
omxCheckCloseEnough(multivariateRegFit$output$estimate[["varx"]], 1.1048, 0.001)
omxCheckCloseEnough(multivariateRegFit$output$estimate[["varz"]], 0.8276, 0.001)
omxCheckCloseEnough(multivariateRegFit$output$estimate[["covxz"]], 0.2861, 0.001)
omxCheckCloseEnough(multivariateRegFit$output$estimate[["meanx"]], 0.0540, 0.001)
omxCheckCloseEnough(multivariateRegFit$output$estimate[["meanz"]], 4.0610, 0.001)

# Compare OpenMx results to Mx results 
# -----------------------------------------------------------------------------
