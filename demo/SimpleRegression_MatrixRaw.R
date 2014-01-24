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


# -----------------------------------------------------------------------------
# Program: SimpleRegression_MatrixRaw.R  
# Author: Ryne Estabrook
# Date: 2009.08.01 
#
# ModelType: Regression
# DataType: Continuous
# Field: None
#
# Purpose:
#      Simple Regression model to estimate effect of 
#      independent on dependent variables
#      Matrix style model input - Raw data input
#
# RevisionHistory:
#      Hermine Maes -- 10 08 2009 updated & reformatted
#      Ross Gore -- 06 06 2011	added Model, Data & Field metadata
# -----------------------------------------------------------------------------

require(OpenMx)
# Load Library
# -----------------------------------------------------------------------------

data(myRegDataRaw)
SimpleDataRaw <- myRegDataRaw[,c("x","y")]
# Prepare Data
# -----------------------------------------------------------------------------

uniRegModel <- mxModel("Simple Regression Matrix Specification", 
    mxData(
        observed=SimpleDataRaw,
        type="raw"
    ),
    mxMatrix(
        type="Full", 
        nrow=2, 
        ncol=2,
        free=c(F, F,
               T, F),
        values=c(0, 0,
                 1, 0),
        labels=c(NA,     NA,
                "beta1", NA),
        byrow=TRUE,
        name="A"
    ),
    mxMatrix(
        type="Symm", 
        nrow=2, 
        ncol=2, 
        values=c(1, 0,
                 0, 1),
        free=c(T, F,
               F, T),
        labels=c("varx", NA,
                  NA,    "residual"),
        byrow=TRUE,
        name="S"
    ),
    mxMatrix(
        type="Iden",  
        nrow=2, 
        ncol=2,
        name="F"
    ),
    mxMatrix(
        type="Full", 
        nrow=1, 
        ncol=2,
        free=c(T, T),
        values=c(0, 0),
        labels=c("meanx", "beta0"),
        name="M"),
    mxRAMObjective("A", "S", "F", "M",
		dimnames=c("x","y"))
)
# Create an MxModel object
# -----------------------------------------------------------------------------
     
uniRegFit<-mxRun(uniRegModel)
# Fit the model with mxRun
# -----------------------------------------------------------------------------

summary(uniRegFit)
# Print a summary of the results
# -----------------------------------------------------------------------------

uniRegFit@output
# Print the outputs populated by mxRun
# -----------------------------------------------------------------------------

omxCheckCloseEnough(uniRegFit@output$estimate[["beta0"]], 2.5478, 0.001)
omxCheckCloseEnough(uniRegFit@output$estimate[["beta1"]], 0.4831, 0.001)
omxCheckCloseEnough(uniRegFit@output$estimate[["residual"]], 0.6652, 0.001)
omxCheckCloseEnough(uniRegFit@output$estimate[["meanx"]], 0.0542, 0.001)
omxCheckCloseEnough(uniRegFit@output$estimate[["varx"]], 1.1053, 0.001)
# Compare OpenMx results to Mx results 
# -----------------------------------------------------------------------------