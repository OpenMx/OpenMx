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
# Program: MultivariateRegression_MatrixRaw.R  
#  Author: Ryne Estabrook
#    Date: 2009.08.01 
#
# ModelType: Regression
# DataType: Continuous
# Field: None
#
# Purpose: 
#      Multivariate Regression model to estimate effect of 
#      independent on dependent variables
#      Matrix style model input - Raw data input
#
# RevisionHistory:
#      Hermine Maes -- 2009.10.08 updated & reformatted
#      Ross Gore -- 2011.06.15 added Model, Data & Field metadata
# -----------------------------------------------------------------------------

require(OpenMx)
# Load Library
# -----------------------------------------------------------------------------

data(myRegDataRaw)
# Prepare Data
# -----------------------------------------------------------------------------

	

multivariateRegModel <- mxModel("Multiple Regression Matrix Specification", 
    mxData(
    	observed=myRegDataRaw,
    	type="raw"
    ),
    mxMatrix(
    	type="Full", 
    	nrow=4, 
    	ncol=4,
        values=c(0,1,0,1,
                 0,0,0,0,
                 0,1,0,1,
                 0,0,0,0),
        free=c(F, T, F, T,
               F, F, F, F,
               F, T, F, T,
               F, F, F, F),
        labels=c(NA, "betawx", NA, "betawz",
                 NA,  NA,     NA,  NA, 
                 NA, "betayx", NA, "betayz",
                 NA,  NA,     NA,  NA),
        byrow=TRUE,
        name="A"
    ),
    mxMatrix(
    	type="Symm", 
    	nrow=4, 
    	ncol=4, 
        values=c(1,  0, 0,  0,
                 0,  1, 0, .5,
                 0,  0, 1,  0,
                 0, .5, 0,  1),
        free=c(T, F, F, F,
               F, T, F, T,
               F, F, T, F,
               F, T, F, T),
        labels=c("residualw",  NA,     NA,         NA,
                  NA,         "varx",  NA,        "covxz",
                  NA,          NA,    "residualy", NA,
                  NA,         "covxz", NA,        "varz"),
        byrow=TRUE,
        name="S"
    ),
    mxMatrix(
    	type="Iden",
    	nrow=4, 
    	ncol=4,
        name="F"
    ),
    mxMatrix(
    	type="Full", 
    	nrow=1, 
    	ncol=4,
        values=c(0,0,0,0),
        free=c(T,T,T,T),
        labels=c("betaw","meanx","betay","meanz"),
        name="M"
    ),
    mxRAMObjective("A","S","F","M",
		dimnames=c("w", "x", "y", "z"))
)
# Create an MxModel object
# -----------------------------------------------------------------------------
      
multivariateRegFit<-mxRun(multivariateRegModel)

summary(multivariateRegFit)
multivariateRegFit@output


omxCheckCloseEnough(multivariateRegFit@output$estimate[["betay"]], 1.6332, 0.001)
omxCheckCloseEnough(multivariateRegFit@output$estimate[["betayx"]], 0.4246, 0.001)
omxCheckCloseEnough(multivariateRegFit@output$estimate[["betayz"]], 0.2260, 0.001)
omxCheckCloseEnough(multivariateRegFit@output$estimate[["residualy"]], 0.6267, 0.001)
omxCheckCloseEnough(multivariateRegFit@output$estimate[["betaw"]], 0.5139, 0.001)
omxCheckCloseEnough(multivariateRegFit@output$estimate[["betawx"]], -0.2310, 0.001)
omxCheckCloseEnough(multivariateRegFit@output$estimate[["betawz"]], 0.5122, 0.001)
omxCheckCloseEnough(multivariateRegFit@output$estimate[["residualw"]], 0.5914, 0.001)
omxCheckCloseEnough(multivariateRegFit@output$estimate[["varx"]], 1.1053, 0.001)
omxCheckCloseEnough(multivariateRegFit@output$estimate[["varz"]], 0.8275, 0.001)
omxCheckCloseEnough(multivariateRegFit@output$estimate[["covxz"]], 0.2862, 0.001)
omxCheckCloseEnough(multivariateRegFit@output$estimate[["meanx"]], 0.0542, 0.001)
omxCheckCloseEnough(multivariateRegFit@output$estimate[["meanz"]], 4.0611, 0.001)
# Compare OpenMx results to Mx results 
# -----------------------------------------------------------------------------