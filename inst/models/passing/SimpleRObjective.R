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

A <- mxMatrix(nrow = 2, ncol = 2, values = c(1:4), free = TRUE, name = 'A')

squared <- function(x) { x ^ 2 }

objFunction <- function(model, state) {
	return(mxEval(squared(A[1,1] - 4) + 
		squared(A[1,2] - 3) +
		squared(A[2,1] - 2) +
		squared(A[2,2] - 1), model))
}
objective <- mxFitFunctionR(objFunction)

model <- mxModel('model', A, objective)

modelOut <- mxRun(model)

omxCheckCloseEnough(mxEval(A, modelOut), 
	rbind(c(4, 3), c(2, 1)), 
	epsilon = 0.001)

omxCheckCloseEnough(modelOut$output$hessian, diag(2, 4), 1e-9)
