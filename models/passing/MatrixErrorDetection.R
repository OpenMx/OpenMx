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
A <- mxMatrix('Full', 1, 1, labels = 'data.foo', free = TRUE, name = 'A')
model <- mxModel('model', A)
omxCheckError(mxRun(model), 
	paste("The definition variable 'data.foo'",
		"has been assigned to a",
		"free parameter in matrix 'A'"))
A <- mxMatrix('Full', 1, 1, labels = 'foo[1,2]', free = TRUE, name = 'A')
model <- mxModel('model', A)
omxCheckError(mxRun(model), 
	paste("The substitution 'foo[1,2]'",
		"has been assigned to a",
		"free parameter in matrix 'A'"))
A <- mxMatrix('Full', 1, 1, labels = 'model2.B[1,1]', name = 'A')
B <- mxMatrix('Full', 1, 1, labels = 'model1.A[1,1]', name = 'B')
model1 <- mxModel('model1', A)
model2 <- mxModel('model2', B)
model <- mxModel('model', model1, model2)
omxCheckError(mxRun(model),
	paste("A cycle has been detected",
		"in model 'model' involving the",
		"following elements: 'model2.B'",
		"and 'model1.A'"))
A <- mxMatrix('Full', 1, 1, labels = 'B[1,1]', name = 'A')
B <- mxMatrix('Full', 1, 1, labels = 'C[1,1]', name = 'B')
C <- mxMatrix('Full', 1, 1, labels = 'A[1,1]', name = 'C')
model <- mxModel('model', A, B, C)
omxCheckError(mxRun(model),
	paste("A cycle has been detected",
		"in model 'model' involving the",
		"following elements: 'B',",
		"'C', and 'A'"))
A <- mxMatrix('Full', 2, 2, name = 'A')
B <- mxMatrix('Full', 2, 2, name = 'B')
C <- mxAlgebra(A + B, name = 'C')
D <- mxMatrix('Full', 1, 1, labels = 'C[2,3]', name = 'D')
model <- mxModel('model', A, B, C, D)
omxCheckError(mxRun(model),
	paste("The substitution 'C[2,3]'",
		"detected in the matrix 'D'",
		"in model 'model' has invalid",
		"(row,col) values"))
