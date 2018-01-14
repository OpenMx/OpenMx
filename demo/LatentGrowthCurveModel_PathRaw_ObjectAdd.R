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
# Program: LatentGrowthModel_PathRaw_ObjectAdd.R  
# Author: Michael Hunter
# Date: 2011.07.22
#
# ModelType: Growth Curve
# DataType: Longitudinal
# Field: None
#
# Purpose: 
#      Latent Growth model to estimate means and 
#      (co)variances of slope and intercept
#      Path style model input - Raw data input
#      Piecemeal
#
# RevisionHistory:
#      Michael Hunter -- 2011.07.22 Took template from Ryne Estabrook's LatentGrowthModel_PathRaw.R
# -----------------------------------------------------------------------------

require(OpenMx)
# Load Library
# -----------------------------------------------------------------------------

data(myLongitudinalData)
# Prepare Data
# -----------------------------------------------------------------------------

# data for model
gcData <- mxData(observed=myLongitudinalData, type="raw")

# residual variances
gcResid <- mxPath(
    	from=c("x1","x2","x3","x4","x5"), 
        arrows=2,
        free=TRUE, 
        values = c(1, 1, 1, 1, 1),
        labels=c("residual","residual","residual","residual","residual")
)

# latent variances and covariance
gcLatCov <- mxPath(
    	from=c("intercept","slope"), 
        arrows=2,
		connect="unique.pairs",
        free=TRUE, 
        values=c(1, 1, 1),
        labels=c("vari", "cov", "vars")
)

# intercept loadings
gcIntercept <- mxPath(
    	from="intercept",
        to=c("x1","x2","x3","x4","x5"),
        arrows=1,
        free=FALSE,
        values=c(1, 1, 1, 1, 1)
)

# slope loadings
gcSlope <- mxPath(
    	from="slope",
        to=c("x1","x2","x3","x4","x5"),
        arrows=1,
        free=FALSE,
        values=c(0, 1, 2, 3, 4)
)
    
# manifest means
gcManMeans <- mxPath(from="one",
        to=c("x1", "x2", "x3", "x4", "x5"),
        arrows=1,
        free=FALSE,
        values=c(0, 0, 0, 0, 0)
)

# latent means
gcLatMeans <- mxPath(from="one",
        to=c("intercept", "slope"),
        arrows=1,
        free=TRUE,
        values=c(1, 1),
        labels=c("meani", "means")
)

growthCurveModel <- mxModel(
    name="Linear Growth Curve Model Path Specification", 
    type="RAM",
    manifestVars=c("x1","x2","x3","x4","x5"),
    latentVars=c("intercept","slope"),
    gcData,
    gcResid,
    gcLatCov,
    gcIntercept,
    gcSlope,
    gcManMeans,
    gcLatMeans
)

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
