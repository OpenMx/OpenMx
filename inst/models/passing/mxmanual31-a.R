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

A <- mxMatrix(values = 0.5, nrow = 2, ncol = 1, 
	free = TRUE, name = "A")

D <- mxMatrix(type = "Diag", values = c(0, 0.5), 
	free = c(FALSE, TRUE), nrow = 2, name = "D")

D$lbound[2,2] <- 0.001
A$lbound[1,1] <- 0.001

expectedCov <- mxAlgebra(A %*% t(A) + D, "expectedCov", 
	dimnames = list(varNames, varNames))

observedCov <- mxData(matrix(c(1.2, 0.8, 0.8, 1.3), 2, 2, 
	dimnames = list(varNames, varNames)), 'cov', numObs = 150)

objective <- mxExpectationNormal(covariance = "expectedCov")

model <- mxModel(A, D, expectedCov, objective, observedCov, mxFitFunctionML())

model <- mxRun(model)

omxCheckCloseEnough(model$output$estimate, c(1.0917, .7278, .7615), 0.001)
