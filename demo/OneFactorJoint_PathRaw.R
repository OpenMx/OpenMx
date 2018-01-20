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
# Program: OneFactorJoint_PathRaw.R  
# Author: Ryne Estabrook
# Date: 2014.05.09 
#
# ModelType: Factor
# DataType: Ordinal
# Field: None
#
# Purpose:
#      One Factor model to estimate factor loadings, residual variances, means and thresholds
#      Path style model input - Raw data input
#
# RevisionHistory:
#      Hermine Maes -- 2014.11.02 piecewise specification
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


#------------------------------------------------------------------------------
# Create an MxModel object
dataRaw      <- mxData( observed=oneFactorJoint, type="raw" )
# residual variances
resVars      <- mxPath( from=c("x1","x2","x3","z1","z2","z3"), arrows=2,
                        free=c(TRUE, TRUE, TRUE, FALSE, FALSE, FALSE),
                        values=1, labels=c("e1","e2","e3","e4","e5","e6") )
# latent variance
latVar       <- mxPath( from="F1", arrows=2,
                        free=FALSE, values=1, labels ="varF1" )
# factor loadings
facLoads     <- mxPath( from="F1", to=c("x1","x2","x3","z1","z2","z3"), arrows=1,
                        free=TRUE, values=1, labels=c("l1","l2","l3","l4","l5","l6") )
# means
means        <- mxPath( from="one", to=c("x1","x2","x3","z1","z2","z3","F1"), arrows=1,
                        free=c(TRUE,TRUE,TRUE,FALSE,FALSE,FALSE,FALSE), values=0,
                        labels=c("meanx1","meanx2","meanx3","meanz1","meanz2","meanz3","meanF") )
# thresholds
thresholds   <- mxThreshold(vars=c("z1","z2","z3"), nThresh=c(1,1,2),
                        free=TRUE, values=c(-1, 0, -.5, 1.2),
                        labels=c("var1_th1","var2_th1","var3_th1","var3_th2") )

oneFactorJointModel <- mxModel("Common Factor Model Path Specification", type="RAM",
                        manifestVars=c("x1","x2","x3","z1","z2","z3"), latentVars="F1",
                        dataRaw, resVars, latVar, means, thresholds, facLoads)

#------------------------------------------------------------------------------
# Fit the model with mxRun
oneFactorJointFit <- mxRun(oneFactorJointModel)

#------------------------------------------------------------------------------
# Print a summary of the results
summary(oneFactorJointFit)
coef(oneFactorJointFit)

#------------------------------------------------------------------------------
# Generate data from the model
simulatedData <- mxGenerateData(oneFactorJointModel, 100)

omxCheckEquals(dim(simulatedData), c(100,6) )

#------------------------------------------------------------------------------
