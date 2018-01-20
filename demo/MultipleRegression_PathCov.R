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
# Program: MultipleRegression_PathCov.R  
# Author: Ryne Estabrook
# Date: 2009.08.01 
#
# ModelType: Regression
# DataType: Continuous
# Field: None
#
# Purpose:
#      Multiple Regression model to estimate effect of 
#      independent on dependent variables
#      Path style model input - Covariance matrix data input
#
# RevisionHistory:
#      Hermine Maes -- 2009.10.08 updated & reformatted
#      Ross Gore -- 2011.06.15 added Model, Data and Field metadata
#      Hermine Maes -- 2014.11.02 piecewise specification
# -----------------------------------------------------------------------------

require(OpenMx)
# Load Library
# -----------------------------------------------------------------------------

myRegDataCov <- matrix(
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

MultipleDataCov <- myRegDataCov[c("x","y","z"),c("x","y","z")]	
MultipleDataMeans <- myRegDataMeans[c(2,3,4)]
# Prepare Data
# -----------------------------------------------------------------------------

# dataset
dataCov      <- mxData( observed=MultipleDataCov,  type="cov", numObs=100, 
                        means=MultipleDataMeans )
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
                        dataCov, manifestVars=c("x","y","z"), 
                        varPaths, covPaths, regPaths, means)
# Create an MxModel object
# -----------------------------------------------------------------------------
      
multiRegFit<-mxRun(multiRegModel)

summary(multiRegFit)
multiRegFit$output


omxCheckCloseEnough(multiRegFit$output$estimate[["beta0"]], 1.6312, 0.001)
omxCheckCloseEnough(multiRegFit$output$estimate[["betax"]], 0.4243, 0.001)
omxCheckCloseEnough(multiRegFit$output$estimate[["betaz"]], 0.2265, 0.001)
omxCheckCloseEnough(multiRegFit$output$estimate[["res"]], 0.6272, 0.001)
omxCheckCloseEnough(multiRegFit$output$estimate[["varx"]], 1.1048, 0.001)
omxCheckCloseEnough(multiRegFit$output$estimate[["varz"]], 0.8276, 0.001)
omxCheckCloseEnough(multiRegFit$output$estimate[["covxz"]], 0.2861, 0.001)
omxCheckCloseEnough(multiRegFit$output$estimate[["meanx"]], 0.0540, 0.001)
omxCheckCloseEnough(multiRegFit$output$estimate[["meanz"]], 4.0610, 0.001)
# Compare OpenMx results to Mx results 
# -----------------------------------------------------------------------------
