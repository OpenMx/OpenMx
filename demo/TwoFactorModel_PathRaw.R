#
#   Copyright 2007-2019 by the individuals mentioned in the source code history
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
# Program: TwoFactorModel_PathRaw.R  
# Author: Ryne Estabrook
# Date: 2009.08.01 
#
# ModelType: Factor
# DataType: Continuous
# Field: None
#
# Purpose: 
# 	   Two Factor model to estimate factor loadings, residual variances and means
# 	   Path style model input - Raw data input
#
# RevisionHistory:
# 	   Hermine Maes -- 2009.10.08 updated & reformatted
# 	   Ross Gore -- 2011.06.06	added Model, Data & Field metadata
#      Hermine Maes -- 2014.11.02 piecewise specification
# -----------------------------------------------------------------------------

require(OpenMx)
# Load Library
# -----------------------------------------------------------------------------

data(myFADataRaw)
# Prepare Data
# -----------------------------------------------------------------------------

twoFactorRaw <- myFADataRaw[,c("x1","x2","x3","y1","y2","y3")]

dataRaw      <- mxData( observed=twoFactorRaw, type="raw" )
# residual variances
resVars      <- mxPath( from=c("x1", "x2", "x3", "y1", "y2", "y3"), arrows=2,
                        free=TRUE, values=c(1,1,1,1,1,1),
                        labels=c("e1","e2","e3","e4","e5","e6") ) 
# latent variances and covariance
latVars      <- mxPath( from=c("F1","F2"), arrows=2, connect="unique.pairs",
                        free=TRUE, values=c(1,.5,1), labels=c("varF1","cov","varF2") )
# factor loadings for x variables	
facLoadsX    <- mxPath( from="F1", to=c("x1","x2","x3"), arrows=1,
                        free=c(F,T,T), values=c(1,1,1), labels=c("l1","l2","l3") )
# factor loadings for y variables
facLoadsY    <- mxPath( from="F2", to=c("y1","y2","y3"), arrows=1,
                        free=c(F,T,T), values=c(1,1,1), labels=c("l4","l5","l6") )
# means
means        <- mxPath( from="one", to=c("x1","x2","x3","y1","y2","y3","F1","F2"), 
                        arrows=1,
                        free=c(T,T,T,T,T,T,F,F), values=c(1,1,1,1,1,1,0,0),
                        labels=c("meanx1","meanx2","meanx3",
                                 "meany1","meany2","meany3",NA,NA) ) 

twoFactorModel <- mxModel("Two Factor Model Path Specification", type="RAM",
                        manifestVars=c("x1", "x2", "x3", "y1", "y2", "y3"), 
                        latentVars=c("F1","F2"),
                        dataRaw, resVars, latVars, facLoadsX, facLoadsY, means)
# Create an MxModel object
# -----------------------------------------------------------------------------
      
twoFactorFit <- mxRun(twoFactorModel)

summary(twoFactorFit)
coef(twoFactorFit)

omxCheckCloseEnough(coef(twoFactorFit)[["l2"]], 0.9723, 0.01)
omxCheckCloseEnough(coef(twoFactorFit)[["l3"]], 0.9313, 0.01)
omxCheckCloseEnough(coef(twoFactorFit)[["l5"]], 1.0498, 0.01)
omxCheckCloseEnough(coef(twoFactorFit)[["l6"]], 1.0531, 0.01)
omxCheckCloseEnough(coef(twoFactorFit)[["varF1"]], 0.6604, 0.01)
omxCheckCloseEnough(coef(twoFactorFit)[["varF2"]], 0.4505, 0.01)
omxCheckCloseEnough(coef(twoFactorFit)[["cov"]], 0.2952, 0.01)
omxCheckCloseEnough(coef(twoFactorFit)[["e1"]], 0.3349, 0.01)
omxCheckCloseEnough(coef(twoFactorFit)[["e2"]], 0.3985, 0.01)
omxCheckCloseEnough(coef(twoFactorFit)[["e3"]], 0.4091, 0.01)
omxCheckCloseEnough(coef(twoFactorFit)[["e4"]], 0.5404, 0.01)
omxCheckCloseEnough(coef(twoFactorFit)[["e5"]], 0.4809, 0.01)
omxCheckCloseEnough(coef(twoFactorFit)[["e6"]], 0.5571, 0.01)
omxCheckCloseEnough(coef(twoFactorFit)[["meanx1"]], 2.988, 0.01)
omxCheckCloseEnough(coef(twoFactorFit)[["meanx2"]], 3.0113, 0.01)
omxCheckCloseEnough(coef(twoFactorFit)[["meanx3"]], 2.9861, 0.01)
omxCheckCloseEnough(coef(twoFactorFit)[["meany1"]], 2.9554, 0.01)
omxCheckCloseEnough(coef(twoFactorFit)[["meany2"]], 2.9562, 0.01)
omxCheckCloseEnough(coef(twoFactorFit)[["meany3"]], 2.9673, 0.01)
# Compare OpenMx results to Mx results 
# -----------------------------------------------------------------------------
