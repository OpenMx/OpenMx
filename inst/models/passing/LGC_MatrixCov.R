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

myLongitudinalDataCov <- matrix(
    c(6.362, 4.344, 4.915,  5.045,  5.966,
  	  4.344, 7.241, 5.825,  6.181,  7.252,
  	  4.915, 5.825, 9.348,  7.727,  8.968,
  	  5.045, 6.181, 7.727, 10.821, 10.135,
  	  5.966, 7.252, 8.968, 10.135, 14.220),
  	nrow=5,
  	dimnames=list(
        c("x1","x2","x3","x4","x5"),
        c("x1","x2","x3","x4","x5"))
	)

myLongitudinalDataMean <- c(9.864, 11.812, 13.612, 15.317, 17.178)
names(myLongitudinalDataMean) <- c("x1","x2","x3","x4","x5")

growthCurveModel <- mxModel("Linear Growth Curve Model, Matrix Specification", 
    mxData(myLongitudinalDataCov, 
        type="cov", 
        numObs=500,
        mean=myLongitudinalDataMean),
    mxMatrix(
        type="Full",
        nrow=7, 
        ncol=7,
        free=F,
        values=c(0,0,0,0,0,1,0,
                 0,0,0,0,0,1,1,
                 0,0,0,0,0,1,2,
                 0,0,0,0,0,1,3,
                 0,0,0,0,0,1,4,
                 0,0,0,0,0,0,0,
                 0,0,0,0,0,0,0),
        byrow=TRUE,
        name="A"),
    mxMatrix(
        type="Symm",
        nrow=7,
        ncol=7,
        free=c(T, F, F, F, F, F, F,
               F, T, F, F, F, F, F,
               F, F, T, F, F, F, F,
               F, F, F, T, F, F, F,
               F, F, F, F, T, F, F,
               F, F, F, F, F, T, T,
               F, F, F, F, F, T, T),
        values=c(0,0,0,0,0,  0,  0,
                 0,0,0,0,0,  0,  0,
                 0,0,0,0,0,  0,  0,
                 0,0,0,0,0,  0,  0,
                 0,0,0,0,0,  0,  0,
                 0,0,0,0,0,  1,0.5,
                 0,0,0,0,0,0.5,  1),
        labels=c("residual", NA, NA, NA, NA, NA, NA,
                 NA, "residual", NA, NA, NA, NA, NA,
                 NA, NA, "residual", NA, NA, NA, NA,
                 NA, NA, NA, "residual", NA, NA, NA,
                 NA, NA, NA, NA, "residual", NA, NA,
                 NA, NA, NA, NA, NA, "vari", "cov",
                 NA, NA, NA, NA, NA, "cov", "vars"),
        byrow= TRUE,
        name="S"),
    mxMatrix(
        type="Full",
        nrow=5,
        ncol=7,
        free=F,
        values=c(1,0,0,0,0,0,0,
                 0,1,0,0,0,0,0,
                 0,0,1,0,0,0,0,
                 0,0,0,1,0,0,0,
                 0,0,0,0,1,0,0),
        byrow=T,
        name="F"),
    mxMatrix("Full",  nrow=1, ncol=7,
        values=c(0,0,0,0,0,1,1),
        free=c(F,F,F,F,F,T,T),
        labels=c(NA,NA,NA,NA,NA,"meani","means"),
        name="M"),
    mxFitFunctionML(),mxExpectationRAM("A","S","F","M",dimnames=c("x1","x2","x3","x4","x5","intercept","slope"))
    )
      
growthCurveFit<-mxRun(growthCurveModel, suppressWarnings=TRUE)

omxCheckCloseEnough(growthCurveFit$output$estimate[["meani"]], 9.930, 0.01)
omxCheckCloseEnough(growthCurveFit$output$estimate[["means"]], 1.813, 0.01)
omxCheckCloseEnough(growthCurveFit$output$estimate[["vari"]], 3.886, 0.01)
omxCheckCloseEnough(growthCurveFit$output$estimate[["vars"]], 0.258, 0.01)
omxCheckCloseEnough(growthCurveFit$output$estimate[["cov"]], 0.460, 0.01)
omxCheckCloseEnough(growthCurveFit$output$estimate[["residual"]], 2.316, 0.01)
