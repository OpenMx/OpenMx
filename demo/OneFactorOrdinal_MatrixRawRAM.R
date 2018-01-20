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
# Program: OneFactorOrdinal_MatrixRawRAM.R  
# Author: Ryne Estabrook
# Date: 2014.05.09 
#
# ModelType: Factor
# DataType: Ordinal
# Field: None
#
# Purpose:
#      One Factor model to estimate factor loadings, residual variances, means and thresholds
#      RAM Matrix style model input - Raw data input
#
# RevisionHistory:
#      Hermine Maes -- 2014.11.02 piecewise specification
# -----------------------------------------------------------------------------

require(OpenMx)
# Load Library
# -----------------------------------------------------------------------------

data(myFADataRaw)

oneFactorOrd <- myFADataRaw[,c("z1","z2","z3")]

oneFactorOrd$z1 <- mxFactor(oneFactorOrd$z1, levels=c(0, 1))
oneFactorOrd$z2 <- mxFactor(oneFactorOrd$z2, levels=c(0, 1))
oneFactorOrd$z3 <- mxFactor(oneFactorOrd$z3, levels=c(0, 1, 2))

# Prepare Data
# -----------------------------------------------------------------------------

dataRaw      <- mxData(oneFactorOrd, type="raw")
# asymmetric paths
matrA        <- mxMatrix( type="Full", nrow=4, ncol=4,
                          free=c(F,F,F,T,
                                 F,F,F,T,
                                 F,F,F,T,
                                 F,F,F,F),
                          values=c(0,0,0,1,
                                   0,0,0,1,
                                   0,0,0,1,
                                   0,0,0,0),
                          labels=c(NA,NA,NA,"l1",
                                   NA,NA,NA,"l2",
                                   NA,NA,NA,"l3",
                                   NA,NA,NA,NA),
                          byrow=TRUE, name="A" )
# symmetric paths
matrS        <- mxMatrix( type="Symm", nrow=4, ncol=4, 
                          free=FALSE, 
                          values=diag(4),
                          labels=c("e1", NA, NA,  NA,
                                    NA,"e2", NA,  NA,
                                    NA,  NA,"e3", NA,
                                    NA,  NA, NA, "varF1"),
                          byrow=TRUE, name="S" )
# filter matrix
matrF        <- mxMatrix( type="Full", nrow=3, ncol=4,
                          free=FALSE, values=c(1,0,0,0,  0,1,0,0,  0,0,1,0),
                          byrow=TRUE, name="F" )
# means
matrM        <- mxMatrix( type="Full", nrow=1, ncol=4,
                          free=FALSE, values=0,
                          labels=c("meanz1","meanz2","meanz3",NA), name="M" )
thresh       <- mxMatrix( type="Full", nrow=2, ncol=3,
                          free=c(TRUE,TRUE,TRUE,FALSE,FALSE,TRUE), 
                          values=c(-1,0,-.5,NA,NA,1.2), byrow=TRUE, name="thresh" )
exp          <- mxExpectationRAM("A","S","F","M", dimnames=c("z1","z2","z3","F1"), 
                          thresholds="thresh", threshnames=c("z1","z2","z3"))
funML        <- mxFitFunctionML()

oneFactorOrdinalModel <- mxModel("Common Factor Model Matrix Specification", 
                                 dataRaw, matrA, matrS, matrF, matrM, thresh, exp, funML)
# Create an MxModel object
# -----------------------------------------------------------------------------
                        
oneFactorOrdinalFit <- mxRun(oneFactorOrdinalModel)
# Fit the model with mxRun
# -----------------------------------------------------------------------------

summary(oneFactorOrdinalFit)
coef(oneFactorOrdinalFit)
# Print a summary of the results
# -----------------------------------------------------------------------------
