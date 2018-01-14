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
# Program: MultipleRegression_MatrixRaw.R  
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
#      Matrix style model input - Raw data input
#
# RevisionHistory:
#      Hermine Maes -- 2009.10.08 updated & reformatted
#      Ross Gore -- 2011.06.15 added Model, Data & Field metadata
#      Hermine Maes -- 2014.11.02 piecewise specification
# -----------------------------------------------------------------------------

require(OpenMx)
# Load Library
# -----------------------------------------------------------------------------

data(myRegDataRaw)

MultipleDataRaw<-myRegDataRaw[,c("x","y","z")]
# Prepare Data
# -----------------------------------------------------------------------------

dataRaw     <- mxData( observed=MultipleDataRaw, type="raw" )
matrA       <- mxMatrix( type="Full", nrow=3, ncol=3,
                         free=  c(F,F,F,  T,F,T,  F,F,F),
                         values=c(0,0,0,  1,0,1,  0,0,0),
                         labels=c(NA,NA,NA, "betax",NA,"betaz", NA,NA,NA),
                         byrow=TRUE, name="A" )
matrS       <- mxMatrix( type="Symm", nrow=3, ncol=3, 
                         free=c(T,F,T,  F,T,F,  T,F,T),
                         values=c(1,0,.5,  0,1,0,  .5,0,1),
                         labels=c("varx",NA,"covxz", NA,"residual",NA, "covxz",NA,"varz"),
                         byrow=TRUE, name="S" )
matrF       <- mxMatrix( type="Iden", nrow=3, ncol=3, name="F" )
matrM       <- mxMatrix( type="Full", nrow=1, ncol=3, 
                         free=c(T,T,T), values=c(0,0,0),
                         labels=c("meanx","beta0","meanz"), name="M" )
exp         <- mxExpectationRAM("A","S","F","M", dimnames=c("x","y","z") )
funML       <- mxFitFunctionML()
multiRegModel <- mxModel("Multiple Regression Matrix Specification", 
                         dataRaw, matrA, matrS, matrF, matrM, exp, funML)
# Create an MxModel object
# -----------------------------------------------------------------------------
      
multiRegFit<-mxRun(multiRegModel)

summary(multiRegFit)
multiRegFit$output


omxCheckCloseEnough(multiRegFit$output$estimate[["beta0"]], 1.6332, 0.001)
omxCheckCloseEnough(multiRegFit$output$estimate[["betax"]], 0.4246, 0.001)
omxCheckCloseEnough(multiRegFit$output$estimate[["betaz"]], 0.2260, 0.001)
omxCheckCloseEnough(multiRegFit$output$estimate[["residual"]], 0.6267, 0.001)
omxCheckCloseEnough(multiRegFit$output$estimate[["varx"]], 1.1053, 0.001)
omxCheckCloseEnough(multiRegFit$output$estimate[["varz"]], 0.8275, 0.001)
omxCheckCloseEnough(multiRegFit$output$estimate[["covxz"]], 0.2862, 0.001)
omxCheckCloseEnough(multiRegFit$output$estimate[["meanx"]], 0.0542, 0.001)
omxCheckCloseEnough(multiRegFit$output$estimate[["meanz"]], 4.0611, 0.001)
# Compare OpenMx results to Mx results 
# -----------------------------------------------------------------------------