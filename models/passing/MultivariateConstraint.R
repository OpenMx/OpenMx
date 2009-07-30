#
#   Copyright 2007-2009 The OpenMx Project
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

data <- mxData(matrix(c(3.6,2.2,0.5,2.2,3.2,2,0.5,2,3.9), nrow = 3),
	type = "cov", numObs = 100)

s <- mxMatrix("Symm", free = TRUE, 
	values = matrix(c(10,9,9,9,10,9,9,9,10), nrow = 3),
	labels = matrix(c("v1", "c12", "c13", 
		"c12", "v2", "c23", 
		"c13", "c23", "v3"), nrow = 3),
	name = "s")
	
c <- mxMatrix("Full", free = FALSE,
	values = matrix(c(0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5), nrow = 3),
	name= "c")
	
bounds <- mxBounds(c("v1", "v2", "v3"), min = 0)

constraint <- mxConstraint("s", ">", "c", name = "constraint")

model <- mxModel("model", data, s, c, 
	constraint, bounds, mxMLObjective("s"))
	
run <- mxRun(model)

omxCheckCloseEnough(mxEvaluate(s, run), run@data@observed, .001)
