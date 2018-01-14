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
W <- mxMatrix("Symm", nrow = 3, ncol = 3, name = "W",
	values = c(1.0, 0.4, 0.3, 0.9, 0.5, 1.1))
H <- mxMatrix("Symm", nrow = 3, ncol = 3, name = "H",
	values = c(1.2, 0.42, 0.3, 1.0, 0.47, 0.9))
M <- mxMatrix("Full", nrow = 3, ncol = 3, name = "M",
	values = c(0.4, 0.05, 0.22, 0.1, 0.3, 0.11, 0.2, 0.12, 0.5))
D <- mxAlgebra(solve(W) %*% M %*% solve(H), name = "D")
model <- mxModel(W, H, M, D)
model <- mxRun(model)
print(model)

omxCheckCloseEnough(model[['D']]$result, 
	rbind(c(0.4302, -0.2854, 0.1497),
		c(-0.3471, 0.7770, -0.5237),
		c(0.1361, -0.4901, 0.7829)), 0.01)


