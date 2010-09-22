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
foo <- mxAlgebra(A + B, 'foo')
A <- mxMatrix('Full', 1, 2, name = 'A')
B <- mxMatrix('Full', 2, 1, name = 'B')
model <- mxModel('model', A, B, foo)
omxCheckError(mxRun(model), paste("The algebra 'foo' in model",
	"'model' generated the error message:",
	"non-conformable arrays"))
cycle <- mxAlgebra(cycle, 'cycle')
model <- mxModel('model', cycle)
omxCheckError(mxRun(model), paste("A cycle has been detected",
	"in model 'model' involving the following elements:",
	"'cycle'"))
foo <- mxAlgebra(bar, 'foo')
bar <- mxAlgebra(foo, 'bar')
model <- mxModel('model', foo, bar)
omxCheckError(mxRun(model), paste("A cycle has been detected",
	"in model 'model' involving the following elements:",
	"'bar' and 'foo'"))
A <- mxMatrix('Full', 1, 1, name = 'A')
B <- mxMatrix('Full', 2, 2, name = 'B')
C <- mxAlgebra(A, 'C')
D <- mxAlgebra(B, 'D')
constraint1 <- mxConstraint(A == B, name = 'constraint1')
constraint2 <- mxConstraint(C == D, name = 'constraint2')
constraint3 <- mxConstraint(1 == B, name = 'constraint3')
model1 <- mxModel('model1', A, B, C, D, constraint1)
model2 <- mxModel('model2', A, B, C, D, constraint2)
model3 <- mxModel('model3', A, B, C, D, constraint3)
omxCheckError(mxRun(model1), paste("The algebras/matrices",
	"'A' and 'B' in model 'model1' are in constraint 'constraint1'",
	"and are not of identical dimensions. The left-hand side is",
	"1 x 1 and the right-hand side is 2 x 2."))
omxCheckError(mxRun(model2), paste("The algebras/matrices",
	"'C' and 'D' in model 'model2' are in constraint 'constraint2'",
	"and are not of identical dimensions. The left-hand side is",
	"1 x 1 and the right-hand side is 2 x 2."))
omxCheckError(mxRun(model3), paste("The algebras/matrices",
	"'1' and 'B' in model 'model3' are in constraint 'constraint3'",
	"and are not of identical dimensions. The left-hand side is",
	"1 x 1 and the right-hand side is 2 x 2."))
A <- mxMatrix('Full', 1, 1, name = 'A')
B <- mxMatrix('Full', 1, 1, name = 'B', labels = 'A[0,0]')
model <- mxModel('model', A, B)
omxCheckError(mxRun(model), paste("The label 'A[0,0]' of matrix 'B'",
	"in model 'model' does not evaluate to a (1 x 1) matrix."))
kevin <- 'bacon'
B <- mxAlgebra(A[kevin, ], name = 'B')
dimnames(A) <- list('Tom', 'Cruise')
model <- mxModel('model', A, B)
omxCheckError(mxRun(model), paste("The matrix 'model.A' does",
	"not contain the row name 'bacon'"))
model <- mxModel('model', mxModel("model2", mxAlgebra(model2.objective, name="Obj"), mxAlgebraObjective("Obj")))
omxCheckError(mxRun(model), paste("A cycle has been detected in model",
	"'model' involving the following elements:",
	"'model2.Obj' and 'model2.objective'"))
