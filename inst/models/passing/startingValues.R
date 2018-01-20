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


library(OpenMx)

covMatrix <- matrix( c(0.77642931, 0.39590663, 0.39590663, 0.49115615), 
	nrow = 2, ncol = 2, byrow = TRUE, dimnames = list(c('a','b'), c('a','b')))

model <- mxModel("missingtest",
	   mxMatrix("Full", values = c(0,0.2,0,0), name="A", nrow=2, ncol=2),
	   mxMatrix("Symm", values = c(0.8,0,0,0.8), name="S", nrow=2, ncol=2, free=TRUE),
	   mxMatrix("Iden", name="F", nrow=2, ncol=2, dimnames = list(c('a','b'), c('a','b'))),
	   mxData(covMatrix, 'cov', numObs = 100),
	   mxExpectationRAM("A", "S", "F"),
	   mxFitFunctionML())

model[["A"]]$free[2,1] <- TRUE
model[["A"]]$values[2,1] <- NA   # oops

omxCheckError(mxRun(model, silent=TRUE),
	      "Parameter 'missingtest.A[2,1]' has no starting value")

model$A$free[,] <- FALSE
model$S$free[,] <- FALSE
model <- omxAssignFirstParameters(model)  # should do nothing with no free parameters


test<-matrix(rnorm(100,0,1),ncol=2)
colnames(test)<-c('v1','v2')
testm<-mxModel(
  type='RAM',
  mxData(test,type='raw'),
  mxMatrix(name='A',ncol=2,nrow=2),
  mxMatrix(name='S',ncol=2,nrow=2,labels=c('freematrix1[1,1]',NA,NA,'freematrix2[1,1]')),
  mxMatrix(name='freematrix1',values=3,free=T,nrow=1,ncol=1,labels='freeparam1'),
  mxMatrix(name='freematrix2',values=2,free=T,nrow=1,ncol=1,labels='freeparam1'),
  mxMatrix(name='M',ncol=2,nrow=1),
  mxMatrix(name='F',type='Diag',,values=1,ncol=2,nrow=2,dimnames=list(c('v1','v2'),c('v1','v2'))),
  mxExpectationRAM(M='M',dimnames=c('v1','v2')),
  mxFitFunctionML()
)
omxCheckEquals(length(omxGetParameters(testm)), 1L)

vn <- omxCheckWarning(mxModel("varName", type="RAM",
                              manifestVars=':-<', latentVars = ":-)"),
                      c("The name ':-)' is illegal because it contains the characters '-' and ':' in FUN(X[[i]], ...)",
                        "The name ':-<' is illegal because it contains the characters '-', ':', and '<' in FUN(X[[i]], ...)"))
