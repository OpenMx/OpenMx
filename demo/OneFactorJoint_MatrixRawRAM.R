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


#------------------------------------------------------------------------------
# Program: OneFactorJoint_MatrixRawRAM.R  
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
#      Michael Hunter -- 2016.10.20 improve format
#------------------------------------------------------------------------------
# Load Library
require(OpenMx)

#------------------------------------------------------------------------------
# Prepare Data
data(myFADataRaw)

oneFactorJoint <- myFADataRaw[,c("x1","x2","x3","z1","z2","z3")]
	
oneFactorJoint$z1 <- mxFactor(oneFactorJoint$z1, levels=c(0, 1))
oneFactorJoint$z2 <- mxFactor(oneFactorJoint$z2, levels=c(0, 1))
oneFactorJoint$z3 <- mxFactor(oneFactorJoint$z3, levels=c(0, 1, 2))

dataRaw <- mxData(observed=oneFactorJoint, type="raw")

#------------------------------------------------------------------------------
# Create an MxModel object

# asymmetric paths
matrA        <- mxMatrix( type="Full", nrow=7, ncol=7,
                          free=c(rep(c(F,F,F,F,F,F,T),6),rep(F,7)),
                          values=c(rep(c(0,0,0,0,0,0,1),6),rep(F,7)),
                          labels=rbind(cbind(matrix(NA,6,6),matrix(paste("l",1:6,sep=""),6,1)),
                           matrix(NA,1,7)),
                          byrow=TRUE, name="A" )
# symmetric paths
labelsS      <- matrix(NA,7,7); diag(labelsS) <- c(paste("e",1:6,sep=""),"varF1")
matrS        <- mxMatrix( type="Symm", nrow=7, ncol=7, 
                          free= rbind(cbind(matrix(as.logical(diag(3)),3,3),matrix(F,3,4)), 
                           matrix(F,4,7)),
                          values=diag(7), labels=labelsS, byrow=TRUE, name="S" )
# filter matrix
matrF        <- mxMatrix( type="Full", nrow=6, ncol=7,
                          free=FALSE, values=cbind(diag(6),matrix(0,6,1)),
                          byrow=TRUE, name="F" )
# means
matrM        <- mxMatrix( type="Full", nrow=1, ncol=7,
                          free=c(T,T,T,F,F,F,F), values=c(1,1,1,0,0,0,0),
                          labels=c("meanx1","meanx2","meanx3","meanz1","meanz2","meanz3",NA),
                          name="M" )
thresh       <- mxMatrix( type="Full", nrow=2, ncol=3,
                          free=c(TRUE,TRUE,TRUE,FALSE,FALSE,TRUE), 
                          values=c(-1,0,-.5,NA,NA,1.2), byrow=TRUE, name="thresh" )
exp          <- mxExpectationRAM("A","S","F","M", 
                                 dimnames=c("x1","x2","x3","z1","z2","z3","F1"), 
                                 thresholds="thresh", threshnames=c("z1","z2","z3"))
funML        <- mxFitFunctionML()

oneFactorJointModel <- mxModel("Common Factor Model Matrix Specification", 
                               dataRaw, matrA, matrS, matrF, matrM, thresh, exp, funML)

#------------------------------------------------------------------------------
# Fit the model with mxRun
oneFactorJointFit <- mxRun(oneFactorJointModel)

#------------------------------------------------------------------------------
# Print a summary of the results
summary(oneFactorJointFit)
coef(oneFactorJointFit)

#------------------------------------------------------------------------------
