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
# Program: SimpleRegression_MatrixCov.R  
# Author: Ryne Estabrook
# Date: 2009.08.01 
#
# ModelType: Regression
# DataType: Continuous
# Field: None
#
# Purpose:
#      Simple Regression model to estimate effect of independent on dependent variables
#      Matrix style model input - Covariance matrix data input
#
# RevisionHistory:
#      Hermine Maes -- 2009.10.08 updated & reformatted
#      Ross Gore -- 2011.06.06	added Model, Data & Field metadata
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
 
SimpleDataCov <- myRegDataCov[c("x","y"),c("x","y")]	 
myRegDataMeans <- c(2.582, 0.054, 2.574, 4.061)
names(myRegDataMeans) <- c("w","x","y","z")

SimpleDataMeans <- myRegDataMeans[c(2,3)]
# Prepare Data
# -----------------------------------------------------------------------------

dataCov      <- mxData( observed=SimpleDataCov, type="cov", numObs=100, 
                          means=SimpleDataMeans )
matrA        <- mxMatrix( type="Full", nrow=2, ncol=2, 
                          free=c(F,F,T,F), values=c(0,0,1,0), 
                          labels=c(NA,NA,"beta1",NA), byrow=TRUE, name="A" )
matrS        <- mxMatrix( type="Symm", nrow=2, ncol=2, 
                          free=c(T,F,F,T), values=c(1,0,0,1), 
                          labels=c("varx",NA,NA,"residual"), byrow=TRUE, name="S" )
matrF        <- mxMatrix( type="Iden", nrow=2, ncol=2, name="F" )
matrM        <- mxMatrix( type="Full", nrow=1, ncol=2, 
                          free=c(T,T), values=c(0,0), 
                          labels=c("meanx","beta0"), name="M")
expRAM       <- mxExpectationRAM("A","S","F","M", dimnames=c("x","y"))
funML        <- mxFitFunctionML()
uniRegModel  <- mxModel("Simple Regression Matrix Specification", 
                        dataCov, matrA, matrS, matrF, matrM, expRAM, funML)
# Create an MxModel object
# -----------------------------------------------------------------------------
      
uniRegFit <- mxRun(uniRegModel)

summary(uniRegFit)
uniRegFit$output


omxCheckCloseEnough(uniRegFit$output$estimate[["beta0"]], 2.54776, 0.001)
omxCheckCloseEnough(uniRegFit$output$estimate[["beta1"]], 0.48312, 0.001)
omxCheckCloseEnough(uniRegFit$output$estimate[["residual"]], 0.672, 0.01)
omxCheckCloseEnough(uniRegFit$output$estimate[["meanx"]], 0.05412, 0.001)
omxCheckCloseEnough(uniRegFit$output$estimate[["varx"]], 1.10484, 0.001)
# Compare OpenMx results to Mx results 
# -----------------------------------------------------------------------------
