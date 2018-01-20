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
# Program: LatentGrowthModel_PathRaw.R  
# Author: Ryne Estabrook
# Date: 2010.09.17 
#
# ModelType: Growth Mixture
# DataType: Continuous
# Field: None
#
# Purpose: 
#      Growth Mixture Model
#      Path style model input - Raw data input
#
# RevisionHistory:
#      Ryne Estabrook -- 2010.10.18 Changed constraints used to identify mixture
#      Ross Gore -- 2011.06.15 added Model, Data & Field metadata
#      Hermine Maes -- 2014.11.02 piecewise specification
# -----------------------------------------------------------------------------

require(OpenMx)
# Load Libraries
# -----------------------------------------------------------------------------


data(myGrowthMixtureData)
#Prepare Data
# -----------------------------------------------------------------------------

matrA        <- mxMatrix( type="Full", nrow=7, ncol=7,
                          free=F, values=rbind(cbind(matrix(0,5,5),
                          matrix(c(rep(1,5),0:4),5,2)),matrix(0,2,7)),
                          byrow=TRUE, name="A" )
labelsS      <- matrix(NA,5,5); diag(labelsS) <- "residual"
matrS        <- mxMatrix( type="Symm", nrow=7, ncol=7,
                          free=rbind(cbind(matrix(as.logical(diag(5)),5,5),
                          matrix(F,5,2)),cbind(matrix(F,2,5),matrix(T,2,2))),
                          values=rbind(cbind(matrix((diag(5)),5,5),
                          matrix(0,5,2)),cbind(matrix(0,2,5),matrix(c(1,.4,.4,1),2,2))),
                          labels=rbind(cbind(labelsS, matrix(NA,5,2)),cbind(matrix(NA,2,5),
                          matrix(c("vari1","cov1","cov1","vars1"),2,2))),
                          byrow= TRUE, name="S" )
matrF        <- mxMatrix( type="Full", nrow=5, ncol=7,
                          free=F, values=cbind(diag(5),matrix(0,5,2)),
                          byrow=T, name="F" )
matrM        <- mxMatrix( type="Full", nrow=1, ncol=7,
                          free=c(F,F,F,F,F,T,T),
                          values=c(0,0,0,0,0,0,-1),
                          labels=c(NA,NA,NA,NA,NA,"meani1","means1"), name="M" )
exp          <- mxExpectationRAM("A","S","F","M",
                          dimnames=c(names(myGrowthMixtureData),"intercept","slope"))
funML        <- mxFitFunctionML(vector=TRUE)
class1       <- mxModel("Class1", matrA, matrS, matrF, matrM, exp, funML)
#Create an MxModel object
# -----------------------------------------------------------------------------

matrS2       <- mxMatrix( type="Symm", nrow=7, ncol=7,
                          free=rbind(cbind(matrix(as.logical(diag(5)),5,5),
                          matrix(F,5,2)),cbind(matrix(F,2,5),matrix(T,2,2))),
                          values=rbind(cbind(matrix((diag(5)),5,5),
                          matrix(0,5,2)),cbind(matrix(0,2,5),matrix(c(1,.5,.5,1),2,2))),
                          labels=rbind(cbind(labelsS, matrix(NA,5,2)),cbind(matrix(NA,2,5),
                          matrix(c("vari2","cov2","cov2","vars2"),2,2))),
                          byrow= TRUE, name="S2" )
matrM2       <- mxMatrix( type="Full", nrow=1, ncol=7,
                          free=c(F,F,F,F,F,T,T),
                          values=c(0,0,0,0,0,5,1),
                          labels=c(NA,NA,NA,NA,NA,"meani2","means2"), name="M2" )
exp          <- mxExpectationRAM("A","S2","F","M2",
                          dimnames=c(names(myGrowthMixtureData),"intercept","slope"))
class2       <- mxModel( class1, name="Class2", matrS2, matrM2, exp )
#Create an MxModel object# -----------------------------------------------------------------------------

# request that individual likelihoods are used
# required for correct parameterization of class probabilities
classP       <- mxMatrix( type="Full", nrow=2, ncol=1, 
                        free=c(TRUE, FALSE), values=1, lbound=0.001, 
                        labels = c("p1","p2"), name="Props" )
# make a matrix of class probabilities
classS       <- mxAlgebra( Props%x%(1/sum(Props)), name="classProbs" )

algFit       <- mxAlgebra( -2*sum(log(classProbs[1,1]%x%Class1.fitfunction 
                           + classProbs[2,1]%x%Class2.fitfunction)), 
                           name="mixtureObj")
fit          <- mxFitFunctionAlgebra("mixtureObj")
dataRaw      <- mxData( observed=myGrowthMixtureData, type="raw" )
      
gmm          <- mxModel("Growth Mixture Model",
                        dataRaw, class1, class2, classP, classS, algFit, fit )     

gmmFit <- mxRun(gmm, suppressWarnings=TRUE)

summary(gmmFit)

gmmFit$classProbs


omxCheckCloseEnough(gmmFit$output$Minus2LogLikelihood, 8739.05, 0.01)
omxCheckCloseEnough(max(mxEval(classProbs, gmmFit)), 0.6009, 0.01)
omxCheckCloseEnough(min(mxEval(classProbs, gmmFit)), 0.3991, 0.01)
# Check results to see if they are within specified bounds
# -----------------------------------------------------------------------------

