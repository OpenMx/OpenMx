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


require(OpenMx)

varNames <- c('x','y')

dataMZ <- mxData(matrix(c(1,.8,.8,1), nrow = 2, ncol=2, 
	dimnames = list(varNames,varNames)), type="cov",
	numObs=100)
dataDZ <- mxData(matrix(c(1,.5,.5,1), nrow = 2, ncol=2,
	dimnames = list(varNames,varNames)), type="cov",
	numObs=100)

X <- mxMatrix("Full",.6,free=TRUE, labels='param1', nrow=1, ncol=1, name="X")
Y <- mxMatrix("Full",.6,free=TRUE, labels='param2', nrow=1, ncol=1, name="Y")
Z <- mxMatrix("Full",.6,free=TRUE, labels='param3', nrow=1, ncol=1, name="Z")
h <- mxMatrix("Full",.5,free=FALSE, nrow=1, ncol=1, name="h")
A <- mxAlgebra(X * t(X), name="A")
C <- mxAlgebra(Y * t(Y), name="C")
E <- mxAlgebra(Z * t(Z), name="E")
cMZ <- mxAlgebra(rbind(cbind(A+C+E,A+C),cbind(A+C,A+C+E)), 
	name="cMZ")
cDZ <- mxAlgebra(rbind(cbind(A+C+E,h%x%A+C),cbind(h%x%A+C,A+C+E)),
	name="cDZ")

objMZ <- mxExpectationNormal("cMZ", dimnames = varNames)
objDZ <- mxExpectationNormal("cDZ", dimnames = varNames)

modelMZ <- mxModel("modelMZ", dataMZ, X,Y,Z,A,C,E,cMZ, objMZ, mxFitFunctionML())
modelDZ <- mxModel("modelDZ", dataDZ, X,Y,Z,h,A,C,E,cDZ, objDZ, mxFitFunctionML())

twin <- mxAlgebra(modelMZ.objective + modelDZ.objective, name="twin")
obj <- mxFitFunctionAlgebra("twin")

model <- mxModel("both", twin, obj, modelMZ, modelDZ)
modelOut <- mxRun(model)

expectedACE <- c(.6, .2, .2) * 99/100
observedACE <- c(modelOut$modelMZ.A$result, 
	modelOut$modelMZ.C$result, modelOut$modelMZ.E$result)

omxCheckCloseEnough(expectedACE, observedACE, epsilon = 10 ^ -4)
