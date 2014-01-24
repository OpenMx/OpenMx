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
# Program: TwoFactorModel_MatrixRaw.R  
# Author: Ryne Estabrook
# Date: 2009.08.01 
#
# ModelType: Factor
# DataType: Continuous
# Field: None
#
# Purpose: 
#      Two Factor model to estimate factor loadings, residual variances and means
#      Matrix style model input - Raw data input
#
# RevisionHistory:
#      Hermine Maes -- 2009.10.08 updated & reformatted
#      Ross Gore -- 2011.06.06	added Model, Data & Field metadata
# -----------------------------------------------------------------------------

require(OpenMx)
# Load Library
# -----------------------------------------------------------------------------

data(myFADataRaw)
# Prepare Data
# -----------------------------------------------------------------------------

manifestVars <- c("x1","x2","x3","y1","y2","y3")
latentVars <- c("F1","F2")

twoFactorRaw <- myFADataRaw[,manifestVars]

twoFactorModel <- mxModel("Two Factor Model Matrix Specification", 
    type="RAM",
    mxData(
    	observed=twoFactorRaw, 
    	type="raw"
    ),
    mxMatrix(
    	type="Full", 
    	nrow=8, 
    	ncol=8,
        values=c(0,0,0,0,0,0,1,0,
                 0,0,0,0,0,0,1,0,
                 0,0,0,0,0,0,1,0,
                 0,0,0,0,0,0,0,1,
                 0,0,0,0,0,0,0,1,
                 0,0,0,0,0,0,0,1,
                 0,0,0,0,0,0,0,0,
                 0,0,0,0,0,0,0,0),
        free=c(F, F, F, F, F, F, F, F,
               F, F, F, F, F, F, T, F,
               F, F, F, F, F, F, T, F,
               F, F, F, F, F, F, F, F,
               F, F, F, F, F, F, F, T,
               F, F, F, F, F, F, F, T,
               F, F, F, F, F, F, F, F,
               F, F, F, F, F, F, F, F),
        labels=c(NA,NA,NA,NA,NA,NA,"l1", NA,
                 NA,NA,NA,NA,NA,NA,"l2", NA,
                 NA,NA,NA,NA,NA,NA,"l3", NA,
                 NA,NA,NA,NA,NA,NA, NA,"l4",
                 NA,NA,NA,NA,NA,NA, NA,"l5",
                 NA,NA,NA,NA,NA,NA, NA,"l6",
                 NA,NA,NA,NA,NA,NA, NA, NA,
                 NA,NA,NA,NA,NA,NA, NA, NA),
        byrow=TRUE,
        name="A"
    ),
    mxMatrix(
    	type="Symm", 
    	nrow=8, 
    	ncol=8, 
        values=c(1,0,0,0,0,0, 0, 0,
                 0,1,0,0,0,0, 0, 0,
                 0,0,1,0,0,0, 0, 0,
                 0,0,0,1,0,0, 0, 0,
                 0,0,0,0,1,0, 0, 0,
                 0,0,0,0,0,1, 0, 0,
                 0,0,0,0,0,0, 1,.5,
                 0,0,0,0,0,0,.5, 1),
        free=c(T, F, F, F, F, F, F, F,
               F, T, F, F, F, F, F, F,
               F, F, T, F, F, F, F, F,
               F, F, F, T, F, F, F, F,
               F, F, F, F, T, F, F, F,
               F, F, F, F, F, T, F, F,
               F, F, F, F, F, F, T, T,
               F, F, F, F, F, F, T, T),
        labels=c("e1", NA,   NA,   NA,   NA,   NA,    NA,    NA,
                 NA, "e2",   NA,   NA,   NA,   NA,    NA,    NA,
                 NA,   NA, "e3",   NA,   NA,   NA,    NA,    NA,
                 NA,   NA,   NA, "e4",   NA,   NA,    NA,    NA,
                 NA,   NA,   NA,   NA, "e5",   NA,    NA,    NA,
                 NA,   NA,   NA,   NA,   NA, "e6",    NA,    NA,
                 NA,   NA,   NA,   NA,   NA,   NA, "varF1", "cov",
                 NA,   NA,   NA,   NA,   NA,   NA, "cov", "varF2"),
        byrow=TRUE,
        name="S"
    ),
    mxMatrix(
    	type="Full",
    	nrow=6, 
    	ncol=8,
        free=F,
        values=c(1,0,0,0,0,0,0,0,
                 0,1,0,0,0,0,0,0,
                 0,0,1,0,0,0,0,0,
                 0,0,0,1,0,0,0,0,
                 0,0,0,0,1,0,0,0,
                 0,0,0,0,0,1,0,0),
        byrow=T,
        name="F"
    ),
    mxMatrix(
    	type="Full", 
    	nrow=1, 
    	ncol=8,
        values=c(1,1,1,1,1,1,0,0),
        free=c(T,T,T,T,T,T,F,F),
        labels=c("meanx1","meanx2","meanx3","meanx4","meanx5","meanx6",NA,NA),
        name="M"
    ),
    mxRAMObjective("A","S","F","M",
		dimnames=c(manifestVars, latentVars))
)
# Create an MxModel object
# -----------------------------------------------------------------------------
      
twoFactorFit <- mxRun(twoFactorModel)

summary(twoFactorFit)
twoFactorFit@output$estimate


omxCheckCloseEnough(twoFactorFit@output$estimate[["l2"]], 0.9723, 0.01)
omxCheckCloseEnough(twoFactorFit@output$estimate[["l3"]], 0.9313, 0.01)
omxCheckCloseEnough(twoFactorFit@output$estimate[["l5"]], 1.0498, 0.01)
omxCheckCloseEnough(twoFactorFit@output$estimate[["l6"]], 1.0531, 0.01)
omxCheckCloseEnough(twoFactorFit@output$estimate[["varF1"]], 0.6604, 0.01)
omxCheckCloseEnough(twoFactorFit@output$estimate[["varF2"]], 0.4505, 0.01)
omxCheckCloseEnough(twoFactorFit@output$estimate[["cov"]], 0.2952, 0.01)
omxCheckCloseEnough(twoFactorFit@output$estimate[["e1"]], 0.3349, 0.01)
omxCheckCloseEnough(twoFactorFit@output$estimate[["e2"]], 0.3985, 0.01)
omxCheckCloseEnough(twoFactorFit@output$estimate[["e3"]], 0.4091, 0.01)
omxCheckCloseEnough(twoFactorFit@output$estimate[["e4"]], 0.5404, 0.01)
omxCheckCloseEnough(twoFactorFit@output$estimate[["e5"]], 0.4809, 0.01)
omxCheckCloseEnough(twoFactorFit@output$estimate[["e6"]], 0.5571, 0.01)
omxCheckCloseEnough(twoFactorFit@output$estimate[["meanx1"]], 2.988, 0.01)
omxCheckCloseEnough(twoFactorFit@output$estimate[["meanx2"]], 3.0113, 0.01)
omxCheckCloseEnough(twoFactorFit@output$estimate[["meanx3"]], 2.9861, 0.01)
omxCheckCloseEnough(twoFactorFit@output$estimate[["meanx4"]], 2.9554, 0.01)
omxCheckCloseEnough(twoFactorFit@output$estimate[["meanx5"]], 2.9562, 0.01)
omxCheckCloseEnough(twoFactorFit@output$estimate[["meanx6"]], 2.9673, 0.01)
# Compare OpenMx results to Mx results 
# -----------------------------------------------------------------------------