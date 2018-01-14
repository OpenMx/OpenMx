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
foo <- mxAlgebra(A + B, 'foo')
A <- mxMatrix('Full', 1, 2, name = 'A')
B <- mxMatrix('Full', 2, 1, name = 'B')
model <- mxModel('model', A, B, foo)
omxCheckError(mxEval(foo, model, compute=TRUE),
	paste("The following error occurred while",
	"evaluating the subexpression 'model.A + model.B'",
	"during the evaluation of 'foo' in model 'model' : non-conformable arrays"))
cycle <- mxAlgebra(cycle, 'cycle')
model <- mxModel('model', cycle)
omxCheckError(mxRun(model),
	      "A cycle has been detected in model 'model' . It involved the following elements: 'cycle'
A common trigger for this error is not providing a name string as the first parameter to mxModel.")
foo <- mxAlgebra(bar, 'foo')
bar <- mxAlgebra(foo, 'bar')
model <- mxModel('model', foo, bar)
omxCheckError(mxRun(model),
	      "A cycle has been detected in model 'model' . It involved the following elements: 'bar' and 'foo'
A common trigger for this error is not providing a name string as the first parameter to mxModel.")
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
omxCheckError(mxRun(model), 
              "Requested improper value (0, 0) from (1, 1) matrix 'model.A'")
kevin <- 'bacon'
B <- mxAlgebra(A[kevin, ], name = 'B')
dimnames(A) <- list('Tom', 'Cruise')
model <- mxModel('model', A, B)
omxCheckError(mxRun(model), paste("The matrix 'model.A' does",
	"not contain the row name 'bacon'"))
model <- mxModel('model', mxModel("model2", mxAlgebra(model2.objective, name="Obj"), mxFitFunctionAlgebra("Obj")))
omxCheckError(mxRun(model),
	      "A cycle has been detected in model 'model' . It involved the following elements: 'model2.Obj' and 'model2.fitfunction'
A common trigger for this error is not providing a name string as the first parameter to mxModel.")
mod <- mxModel("amodel", mxMatrix("Full", 4, 1, values=7, name="M"), mxMatrix("Full", 4, 1, values=1:4, name="Thr"))
omxCheckError(mxAlgebra(expression="minG", name="blah"), "mxAlgebra wants an unquoted expression or formula")
omxCheckWarning(omxMnor(matrix(c(1,90,90,1),2,2), c(0, 0), c(-Inf, -Inf), c(1.282,1.282)), "Correlation with absolute value greater than one found.")

