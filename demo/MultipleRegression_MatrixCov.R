#
#   Copyright 2007-2012 The OpenMx Project
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
# Program: MultipleRegression_MatrixCov.R  
# Author: Ryne Estabrook
# Date: 2009.08.01 
#
# ModelType: Regression
# DataType: Continuous
# Field: None
#
# Purpose: 
#      Multiple Regression model to estimate effect of independent 
#      on dependent variables
#      Matrix style model input - Covariance matrix data input
#
# RevisionHistory:
#      Hermine Maes -- 2009.10.08 updated & reformatted
#      Ross Gore -- 2011.06.15 added Model, Data & Field
# -----------------------------------------------------------------------------

require(OpenMx)
# Load Libraries
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

multiRegModel<-mxModel("Multiple Regression Matrix Specification", 
      mxData(
      		observed=MultipleDataCov, 
      		type="cov", 
      		numObs=100, 
      		mean=MultipleDataMeans
      ),
      mxMatrix("Full", nrow=3, ncol=3,
            values=c(0,0,0,
                     1,0,1,
                     0,0,0),
            free=c(F, F, F,
                   T, F, T,
                   F, F, F),
            labels=c(NA,     NA, NA,
                    "betax", NA,"betaz",
                     NA,     NA, NA),
            byrow=TRUE,
            name="A"
      ),
      mxMatrix("Symm", nrow=3, ncol=3, 
            values=c(1, 0, .5,
                     0, 1, 0,
                    .5, 0, 1),
            free=c(T, F, T,
                   F, T, F,
                   T, F, T),
            labels=c("varx", NA,         "covxz",
                     NA,     "residual",   NA,
                     "covxz",   NA,         "varz"),
            byrow=TRUE,
            name="S"
      ),
      mxMatrix(
      		type="Iden", 
      		nrow=3, 
      		ncol=3,
            name="F"
      ),
      mxMatrix(
      		type="Full", 
      		nrow=1, 
      		ncol=3,
            values=c(0,0,0),
            free=c(T,T,T),
            labels=c("meanx","beta0","meanz"),
            name="M"
      ),
      mxRAMObjective("A","S","F","M",dimnames=c('x','y','z'))
)
# Create an MxModel object
# -----------------------------------------------------------------------------
      
multiRegFit <- mxRun(multiRegModel)

summary(multiRegFit)
multiRegFit@output


omxCheckCloseEnough(multiRegFit@output$estimate[["beta0"]], 1.6312, 0.001)
omxCheckCloseEnough(multiRegFit@output$estimate[["betax"]], 0.4243, 0.001)
omxCheckCloseEnough(multiRegFit@output$estimate[["betaz"]], 0.2265, 0.001)
omxCheckCloseEnough(multiRegFit@output$estimate[["residual"]], 0.6336, 0.001)
omxCheckCloseEnough(multiRegFit@output$estimate[["varx"]], 1.1160, 0.001)
omxCheckCloseEnough(multiRegFit@output$estimate[["varz"]], 0.8360, 0.001)
omxCheckCloseEnough(multiRegFit@output$estimate[["covxz"]], 0.2890, 0.001)
omxCheckCloseEnough(multiRegFit@output$estimate[["meanx"]], 0.0540, 0.001)
omxCheckCloseEnough(multiRegFit@output$estimate[["meanz"]], 4.0610, 0.001)
# Compare OpenMx results to Mx results 
# -----------------------------------------------------------------------------