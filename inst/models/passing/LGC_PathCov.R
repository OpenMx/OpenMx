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

myLongitudinalDataCov<-matrix(
	c(6.362, 4.344, 4.915,  5.045,  5.966,
	  4.344, 7.241, 5.825,  6.181,  7.252,
	  4.915, 5.825, 9.348,  7.727,  8.968,
	  5.045, 6.181, 7.727, 10.821, 10.135,
	  5.966, 7.252, 8.968, 10.135, 14.220),
	nrow=5,
	dimnames=list(
		c("x1","x2","x3","x4","x5"),c("x1","x2","x3","x4","x5"))
	)

myLongitudinalDataMean <- c(9.864, 11.812, 13.612, 15.317, 17.178)
names(myLongitudinalDataMean) <- c("x1","x2","x3","x4","x5")

growthCurveModel <- mxModel("Linear Growth Curve Model, Path Specification", 
    type="RAM",
    mxData(myLongitudinalDataCov, 
        type="cov", 
        numObs=500,
        means=myLongitudinalDataMean),
    manifestVars=c("x1","x2","x3","x4","x5"),
    latentVars=c("intercept","slope"),
    # residual variances
    mxPath(from=c("x1","x2","x3","x4","x5"), 
        arrows=2,
        free=TRUE, 
        values = c(1, 1, 1, 1, 1),
        labels=c("residual","residual","residual","residual","residual")),
    # latent variances and covariance
    mxPath(from=c("intercept","slope"), 
        arrows=2,
		connect="unique.pairs",
        free=TRUE,
		values=c(1, 1, 1),
		labels=c("vari", "cov", "vars")), 
    # intercept loadings
    mxPath(from="intercept",
        to=c("x1","x2","x3","x4","x5"),
        arrows=1,
        free=FALSE,
        values=c(1, 1, 1, 1, 1)),
    # slope loadings
    mxPath(from="slope",
        to=c("x1","x2","x3","x4","x5"),
        arrows=1,
        free=FALSE,
        values=c(0, 1, 2, 3, 4)),
    # manifest means
    mxPath(from="one",
        to=c("x1", "x2", "x3", "x4", "x5"),
        arrows=1,
        free=FALSE,
        values=c(0, 0, 0, 0, 0)),
    # latent means
    mxPath(from="one",
        to=c("intercept", "slope"),
        arrows=1,
        free=TRUE,
        values=c(1, 1),
        labels=c("meani", "means"))
    ) # close model
      
growthCurveFit<-mxRun(growthCurveModel)

omxCheckCloseEnough(growthCurveFit$output$estimate[["meani"]], 9.930, 0.01)
omxCheckCloseEnough(growthCurveFit$output$estimate[["means"]], 1.813, 0.01)
omxCheckCloseEnough(growthCurveFit$output$estimate[["vari"]], 3.886, 0.01)
omxCheckCloseEnough(growthCurveFit$output$estimate[["vars"]], 0.258, 0.01)
omxCheckCloseEnough(growthCurveFit$output$estimate[["cov"]], 0.460, 0.01)
omxCheckCloseEnough(growthCurveFit$output$estimate[["residual"]], 2.316, 0.01)
