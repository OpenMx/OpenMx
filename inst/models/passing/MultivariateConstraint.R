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

varNames <- c('x','y','z')

origCov <- matrix(c(3.6,2.2,0.5,2.2,3.2,2,0.5,2,3.9), nrow = 3,
	dimnames = list(varNames, varNames))

data <- mxData(origCov, type = "cov", numObs = 10)

s <- mxMatrix("Symm", free = TRUE, 
	values = matrix(c(10,9,9,9,10,9,9,9,10), nrow = 3),
	labels = matrix(c("v1", "c12", "c13", 
		"c12", "v2", "c23", 
		"c13", "c23", "v3"), nrow = 3),
	dimnames = list(varNames, varNames),
	name = "s", lbound = -100, ubound = 100)
	
c <- mxMatrix("Full", free = FALSE, values = 0.4, nrow = 3, ncol=3, name= "c")
	
model <- mxModel("model", data, s, c, 
		 mxConstraint(s > c, name = "constraint"),
		 mxBounds(c("v1", "v2", "v3"), min = 0),
		 mxFitFunctionML(),
		 mxExpectationNormal("s"))
	
run <- mxRun(model)

omxCheckCloseEnough((mxEval(s, run)/ origCov), matrix(.9,3,3), 1e-4) # constraint not active

c <- mxMatrix("Full", free = FALSE, values = 0.6, nrow = 3, ncol=3, name= "c")
	
model <- mxModel("model", data, s, c, 
		 mxConstraint(c < s, name = "constraint"),
		 mxBounds(c("v1", "v2", "v3"), min = 0),
		 mxFitFunctionML(),
		 mxExpectationNormal("s"))
	
run <- mxRun(model)

# constraint active
omxCheckCloseEnough(var(c(mxEval(s, run) / origCov)), 0.01381, 1e-4)
#cat(deparse(round(c(vech(mxEval(s, run))),3)))
omxCheckCloseEnough(vech(mxEval(s, run)), c(3.289, 2.073, 0.6, 2.978, 1.906, 3.563), .001)
