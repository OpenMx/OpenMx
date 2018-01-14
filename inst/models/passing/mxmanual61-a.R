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
A <- mxMatrix("Symm", 
	values = c(1.0, 0.23, 0.34, 2.0, 0.45, 3.0),
	nrow = 3, ncol = 3, name = "A")
B <- mxAlgebra(solve(A), name = "B")
model <- mxModel(A, B)
model <- mxRun(model)
print(model)

omxCheckCloseEnough(model[['B']]$result, 
	rbind(c(1.0583, -0.097, -0.1052),
		c(-0.097, 0.5265, -0.0679),
		c(-0.1052, -0.0679, .3554)), 0.01)