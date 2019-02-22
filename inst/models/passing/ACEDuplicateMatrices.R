#
#   Copyright 2007-2018 by the individuals mentioned in the source code history
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

# ====================================================================================================
# = Make an ACE model with duplicated matrices in the MZ and DZ models instead of a shared top model =
# ====================================================================================================

# 1. Make some twin cov data: .8 correlation in MZ, .5 in DZ

varNames <- c('x_T1','x_T2')

dataMZ <- mxData(matrix(c(1,.8,.8,1), nrow = 2, ncol=2, 
	dimnames = list(varNames,varNames)), type="cov", numObs=100)
dataDZ <- mxData(matrix(c(1,.5,.5,1), nrow = 2, ncol=2,
	dimnames = list(varNames,varNames)), type="cov", numObs=100)

# 2. Make a univariate ACE model
h = mxMatrix(name="h", "Full", .5, free=FALSE, nrow=1, ncol=1)
a = mxMatrix(name="a", "Full", .6, free=TRUE, labels='a_r1c1', nrow=1, ncol=1)
c = mxMatrix(name="c", "Full", .6, free=TRUE, labels='c_r1c1', nrow=1, ncol=1)
e = mxMatrix(name="e", "Full", .6, free=TRUE, labels='e_r1c1', nrow=1, ncol=1)
A = mxAlgebra(a * t(a), name="A")
C = mxAlgebra(c * t(c), name="C")
E = mxAlgebra(e * t(e), name="E")
cMZ = mxAlgebra(name="cMZ", rbind(cbind(A+C+E,A+C)    ,cbind(A+C,A+C+E)))
cDZ = mxAlgebra(name="cDZ", rbind(cbind(A+C+E,h%x%A+C),cbind(h%x%A+C,A+C+E)))

objMZ <- mxExpectationNormal("cMZ", dimnames = varNames)
objDZ <- mxExpectationNormal("cDZ", dimnames = varNames)

# this is not great style: no need to duplicate the matrices in each group
MZ <- mxModel("MZ", dataMZ, a,c,e, A,C,E, cMZ  , objMZ, mxFitFunctionML())
DZ <- mxModel("DZ", dataDZ, a,c,e, A,C,E, cDZ,h, objDZ, mxFitFunctionML())

model <- mxModel("both", MZ, DZ, mxFitFunctionMultigroup(c("MZ", "DZ")))
m1 <- mxRun(model)
summary(m1)$parameters

# 3. Derive expectations for A, C, and E, based on .8 and .5 correlations in
# MZ and DZ groups, and check we met them
# A = 2* (.8-.5) = .6
# C = 1 - A+E    = .2
# E = 1 - .8     = .2
# TODO: why * .99?
expectedACE <- c(.6, .2, .2) * 99/100
observedACE <- c(m1$MZ.A$result, m1$MZ.C$result, m1$MZ.E$result)

omxCheckCloseEnough(expectedACE, observedACE, epsilon = 10 ^ -4)
