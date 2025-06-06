#
#   Copyright 2007-2019 by the individuals mentioned in the source code history
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

# -----------------------------------------------------------------------
# Program: LatentGrowthModel_MatrixRaw.R  
# Author: Ryne Estabrook
# Date: 2009.08.01 
#
# ModelType: Growth Curve
# DataType: Longitudinal
# Field: None
#
# Purpose: 
#      Latent Growth model to estimate means and (co)variances of slope and intercept
#      Matrix style model input - Raw data input
#
# RevisionHistory:
#      Hermine Maes -- 2009.10.08 updated & reformatted
#      Ross Gore -- 2011.06.15 added Model, Data & Field metadata 
#      Hermine Maes -- 2014.11.02 piecewise specification
# -----------------------------------------------------------------------

require(OpenMx)
# Load Libraries
# -----------------------------------------------------------------------

data(myLongitudinalData)
# Prepare Data
# -----------------------------------------------------------------------

dataRaw      <- mxData( observed=myLongitudinalData, type="raw" )
matrA        <- mxMatrix( type="Full", nrow=7, ncol=7,
                          free = FALSE,
                          values=c(0,0,0,0,0,1,0,
                                   0,0,0,0,0,1,1,
                                   0,0,0,0,0,1,2,
                                   0,0,0,0,0,1,3,
                                   0,0,0,0,0,1,4,
                                   0,0,0,0,0,0,0,
                                   0,0,0,0,0,0,0),
                          byrow=TRUE, name="A" )
matrS        <- mxMatrix( type="Symm", nrow=7, ncol=7,
                          free=  c(T,F,F,F,F,F,F,
                                   F,T,F,F,F,F,F,
                                   F,F,T,F,F,F,F,
                                   F,F,F,T,F,F,F,
                                   F,F,F,F,T,F,F,
                                   F,F,F,F,F,T,T,
                                   F,F,F,F,F,T,T),
                          values=c(0,0,0,0,0, 0, 0,
                                   0,0,0,0,0, 0, 0,
                                   0,0,0,0,0, 0, 0,
                                   0,0,0,0,0, 0, 0,
                                   0,0,0,0,0, 0, 0,
                                   0,0,0,0,0, 1,.5,
                                   0,0,0,0,0,.5, 1),
                          labels=c("residual", NA, NA, NA, NA, NA, NA,
                                   NA, "residual", NA, NA, NA, NA, NA,
                                   NA, NA, "residual", NA, NA, NA, NA,
                                   NA, NA, NA, "residual", NA, NA, NA,
                                   NA, NA, NA, NA, "residual", NA, NA,
                                   NA, NA, NA, NA, NA, "vari", "cov",
                                   NA, NA, NA, NA, NA, "cov", "vars"),
                          byrow= TRUE, name="S" )
matrF        <- mxMatrix( type="Full", nrow=5, ncol=7,
                          free = FALSE,
                          values=c(1,0,0,0,0,0,0,
                                   0,1,0,0,0,0,0,
                                   0,0,1,0,0,0,0,
                                   0,0,0,1,0,0,0,
                                   0,0,0,0,1,0,0),
                          byrow=TRUE, name="F" )
matrM        <- mxMatrix( type="Full", nrow=1, ncol=7,
                          free=c(F,F,F,F,F,T,T), values=c(0,0,0,0,0,1,1),
                          labels=c(NA,NA,NA,NA,NA,"meani","means"), name="M" )
exp          <- mxExpectationRAM("A","S","F","M", 
                                  dimnames=c(names(myLongitudinalData),"intercept","slope"))
funML        <- mxFitFunctionML()
growthCurveModel <- mxModel("Linear Growth Curve Model Matrix Specification", 
                            dataRaw, matrA, matrS, matrF, matrM, exp, funML)
# Create an MxModel object
# -----------------------------------------------------------------------
      
growthCurveFit<-mxRun(growthCurveModel, suppressWarnings=TRUE)

summary(growthCurveFit)
coef(growthCurveFit)


omxCheckCloseEnough(coef(growthCurveFit)[["meani"]], 9.930, 0.01)
omxCheckCloseEnough(coef(growthCurveFit)[["means"]], 1.813, 0.01)
omxCheckCloseEnough(coef(growthCurveFit)[["vari"]], 3.886, 0.01)
omxCheckCloseEnough(coef(growthCurveFit)[["vars"]], 0.258, 0.01)
omxCheckCloseEnough(coef(growthCurveFit)[["cov"]], 0.460, 0.01)
omxCheckCloseEnough(coef(growthCurveFit)[["residual"]], 2.316, 0.01)
# Compare OpenMx results to Mx results 
# -----------------------------------------------------------------------

# Test `mxFactorScores()`:
fsml <- mxFactorScores(growthCurveFit)
head(fsml)
if(0){
	# TODO why this result??:
	fswml <- mxFactorScores(growthCurveFit,"WeightedML")
	head(fswml)
}
fsr <- mxFactorScores(growthCurveFit,"Regression")
head(fsr)
