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

mxStandardizeRAMpaths(growthCurveFit,T)
rg <- imxRowGradients(growthCurveFit)
rcov <- OpenMx::"%&%"(solve(growthCurveFit$output$hessian/2), nrow(rg)*var(rg/-2))
mxStandardizeRAMpaths(growthCurveFit,SE=T,cov=rcov)


#Check warnings and errors:
pointlessConstraint <- mxModel(growthCurveFit, mxConstraint(1==1))
pointlessConstraint <- mxRun(pointlessConstraint)
omxCheckWarning(
	mxStandardizeRAMpaths(pointlessConstraint, T, rcov),
	"standard errors may be invalid because model 'Linear Growth Curve Model Path Specification' contains at least one mxConstraint"
)
anAlg <- mxAlgebra(1+1)
omxCheckError(
	mxStandardizeRAMpaths(growthCurveFit,T,anAlg),
	"non-NULL value to argument 'cov' must be (or be coercible to) a matrix"
)
omxCheckError(
	mxStandardizeRAMpaths(growthCurveFit,T,matrix(1:6,nrow=2)),
	"non-NULL value to argument 'cov' must be a square matrix; it has 2 rows and 3 columns"
)
rcov2 <- rcov
rownames(rcov2) <- NULL
omxCheckError(
	mxStandardizeRAMpaths(growthCurveFit, T, rcov2),
	"non-NULL value to argument 'cov' must have matching and complete rownames and colnames"
)
omxCheckError(
	mxStandardizeRAMpaths(growthCurveFit, T, rcov[-1,-1]),
	"value of argument 'cov' has dimension 5, but 'Linear Growth Curve Model Path Specification' has 6 free parameters"
)
rcov3 <- rcov
colnames(rcov3) <- c("res","vari","cov","vars","meani","means")
rownames(rcov3) <- c("res","vari","cov","vars","meani","means")
omxCheckError(
	mxStandardizeRAMpaths(growthCurveFit, T, rcov3),
	"the dimnames of the matrix provided for argument 'cov' do not match the free-parameter labels of 'Linear Growth Curve Model Path Specification'"
)
