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
# Program: LatentGrowthModel_PathRaw.R  
# Author: Ryne Estabrook
# Date: 2009.08.01 
#
# ModelType: Growth Curve
# DataType: Longitudinal
# Field: None
#
# Purpose: 
#      Latent Growth model to estimate means and 
#      (co)variances of slope and intercept
#      Path style model input - Raw data input
#
# RevisionHistory:
#      Hermine Maes -- 2009.10.08 updated & reformatted
#      Ross Gore -- 2011.06.15 added Model, Data & Field
#      Hermine Maes -- 2014.11.02 piecewise specification
# -----------------------------------------------------------------------------

require(OpenMx)
# Load Library
# -----------------------------------------------------------------------------

data(myLongitudinalData)
# Prepare Data
# -----------------------------------------------------------------------------

dataRaw      <- mxData( observed=myLongitudinalData, type="raw" )
# residual variances
resVars      <- mxPath( from=c("x1","x2","x3","x4","x5"), arrows=2,
                        free=TRUE,  values = c(1,1,1,1,1),
                        labels=c("residual","residual","residual","residual","residual") )
# latent variances and covariance
latVars      <- mxPath( from=c("intercept","slope"), arrows=2, connect="unique.pairs",
                        free=TRUE, values=c(1,1,1), labels=c("vari","cov","vars") )
# intercept loadings
intLoads     <- mxPath( from="intercept", to=c("x1","x2","x3","x4","x5"), arrows=1,
                        free=FALSE, values=c(1,1,1,1,1) )
# slope loadings
sloLoads     <- mxPath( from="slope", to=c("x1","x2","x3","x4","x5"), arrows=1,
                        free=FALSE, values=c(0,1,2,3,4) )
# manifest means
manMeans     <- mxPath( from="one", to=c("x1","x2","x3","x4","x5"), arrows=1,
                        free=FALSE, values=c(0,0,0,0,0) )
# latent means 
latMeans     <- mxPath( from="one", to=c("intercept", "slope"), arrows=1,
                        free=TRUE, values=c(1,1), labels=c("meani","means") )
growthCurveModel <- mxModel("Linear Growth Curve Model Path Specification", 
                             type="RAM",
                             manifestVars=c("x1","x2","x3","x4","x5"),
                             latentVars=c("intercept","slope"),
                             dataRaw, resVars, latVars, intLoads, sloLoads, 
                             manMeans, latMeans) 
#Create an MxModel object
# -----------------------------------------------------------------------------
      
growthCurveFit <- mxRun(growthCurveModel, suppressWarnings=TRUE)

summary(growthCurveFit)
coef(growthCurveFit)


omxCheckCloseEnough(growthCurveFit$output$estimate[["meani"]], 9.930, 0.01)
omxCheckCloseEnough(growthCurveFit$output$estimate[["means"]], 1.813, 0.01)
omxCheckCloseEnough(growthCurveFit$output$estimate[["vari"]], 3.886, 0.01)
omxCheckCloseEnough(growthCurveFit$output$estimate[["vars"]], 0.258, 0.01)
omxCheckCloseEnough(growthCurveFit$output$estimate[["cov"]], 0.460, 0.01)
omxCheckCloseEnough(growthCurveFit$output$estimate[["residual"]], 2.316, 0.01)
# Compare OpenMx results to Mx results 
# -----------------------------------------------------------------------------
