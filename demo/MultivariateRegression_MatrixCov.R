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
# Program: MultivariateRegression_MatrixCov.R  
# Author: Ryne Estabrook
# Date: 2009.08.01 
#
# ModelType: Regression
# DataType: Continuous
# Field: None
#
# Purpose:
#      Multivariate Regression model to estimate effect of independent on dependent variables
#      Matrix style model input - Covariance matrix data input
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

dataCov     <- mxData( observed=myRegDataCov, type="cov", numObs=100,
		                 mean=myRegDataMeans )
matrA       <- mxMatrix( type="Full", nrow=4, ncol=4,
                         free=  c(F,T,F,T,  F,F,F,F,  F,T,F,T,  F,F,F,F),
                         values=c(0,1,0,1,  0,0,0,0,  0,1,0,1,  0,0,0,0),
                         labels=c(NA,"betawx",NA,"betawz",
                                  NA, NA,     NA, NA, 
                                  NA,"betayx",NA,"betayz",
                                  NA, NA,     NA, NA),
                         byrow=TRUE, name="A" )
matrS       <- mxMatrix( type="Symm", nrow=4, ncol=4, 
                         free=c(T,F,F,F,  F,T,F,T,  F,F,T,F,  F,T,F,T),
                         values=c(1, 0,0, 0,  0, 1,0,.5,  0, 0,1, 0,  0,.5,0, 1),
                         labels=c("residualw", NA,     NA,         NA,
                                   NA,        "varx",  NA,        "covxz",
                                   NA,         NA,    "residualy", NA,
                                   NA,        "covxz", NA,        "varz"),
                         byrow=TRUE, name="S" )
matrF       <- mxMatrix( type="Iden", nrow=4, ncol=4, name="F" )
matrM       <- mxMatrix( type="Full", nrow=1, ncol=4, 
                         free=c(T,T,T,T), values=c(0,0,0,0),
                         labels=c("betaw","meanx","betay","meanz"), name="M" )
exp         <- mxExpectationRAM("A","S","F","M", dimnames=c("w","x","y","z") )
funML       <- mxFitFunctionML()
multivariateRegModel <- mxModel("Multiple Regression Matrix Specification", 
                         dataCov, matrA, matrS, matrF, matrM, exp, funML)
# Create an MxModel object
# -----------------------------------------------------------------------------
      
multivariateRegFit<-mxRun(multivariateRegModel)

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
