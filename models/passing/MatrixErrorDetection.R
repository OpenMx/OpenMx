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
	paste("The label with square brackets",
		"has been assigned to a",
		"free parameter in matrix 'A'",
		"at row 1 and column 1"))
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
	paste("The following error occurred while",
	"evaluating the subexpression 'model.C[2, 3]' during",
	"the evaluation of 'label at row 1 and column 1 of matrix 'D''",
	"in model 'model' : subscript out of bounds"))
# Error check for values parameter
loadings <- c(1, -0.625, 0.1953125, 1,  "h", 0.1953125)
loadings <- matrix(loadings, 2, 3, byrow = TRUE)
omxCheckError(mxMatrix("Full", free=FALSE, values=loadings, name="L", byrow=TRUE),
        paste("'values' argument to mxMatrix function must be of numeric type",
              "in mxMatrix(\"Full\", free = FALSE, values = loadings,",
              "name = \"L\", byrow = TRUE)"))
# Error check for free parameter
loadings <- c(1, -0.625, 0.1953125, 1,  -0.625, 0.1953125)
loadings <- matrix(loadings, 2, 3, byrow = TRUE)
omxCheckError(mxMatrix("Full", free="FALSE", values=loadings, name="L", byrow=TRUE),
        paste("'free' argument to mxMatrix function must be of logical type",
              "in mxMatrix(\"Full\", free = \"FALSE\", values = loadings,",
              "name = \"L\", byrow = TRUE)"))
# Error check for type parameter
loadings <- c(1, -0.625, 0.1953125, 1,  -0.625, 0.1953125)
loadings <- matrix(loadings, 2, 3, byrow = TRUE)
omxCheckError(mxMatrix("wrong", free=FALSE, values=loadings, name="L", byrow=TRUE),
        paste("Type must be one of: Diag Full Iden Lower Stand Sdiag Symm Unit Zero",
              "in mxMatrix(\"wrong\", free = FALSE, values = loadings,",
              "name = \"L\", byrow = TRUE)"))
#Error check both nrow and ncol are na
loadings <- c(1, -0.625, 0.1953125, 1,  -0.625, 0.1953125)
omxCheckError(mxMatrix("Full", free=FALSE, values=loadings, name="L", byrow=TRUE),
        paste("You must specify 'nrow' and 'ncol' arguments in mxMatrix(\"Full\",",
              "free = FALSE, values = loadings, name = \"L\", byrow = TRUE)"))
#Error check upper values of a lower matrix
loadings <- c(1, -0.625, 0.1953125, 1,  -0.625, 0.1953125, 1,  -0.625, 0.1953125)
loadings <- matrix(loadings, 3, 3, byrow = TRUE)
omxCheckError(mxMatrix("Lower", free=FALSE, values=loadings, name="L", byrow=TRUE),
        paste("Upper triangle of values matrix in lower matrix 'L' is not all zeros!",
              "mxMatrix(\"Lower\", free = FALSE, values = loadings,", 
              "name = \"L\", byrow = TRUE)"))
#Error check squareness of a Lower matrix
loadings <- c(1, -0.625, 0.1953125, 1,  -0.625, 0.1953125)
loadings <- matrix(loadings, 2, 3, byrow = TRUE)
omxCheckError(mxMatrix("Lower", free=FALSE, values=loadings, name="L", byrow=TRUE),
		paste("Non-square matrix attempted for lower matrix constructor",
			  "mxMatrix(\"Lower\", free = FALSE, values = loadings,", 
			  "name = \"L\", byrow = TRUE)"))
#Error check squareness of a Symmetrical matrix
loadings <- c(1, -0.625, 0.1953125, 1,  -0.625, 0.1953125)
loadings <- matrix(loadings, 2, 3, byrow = TRUE)
omxCheckError(mxMatrix("Symm", free=FALSE, values=loadings, name="L", byrow=TRUE),
	    paste("Non-square matrix attempted for symmmetric matrix constructor",
			  "mxMatrix(\"Symm\", free = FALSE, values = loadings,", 
			  "name = \"L\", byrow = TRUE)"))
#Error check squareness of a Diagonal matrix
loadings <- c(1, -0.625, 0.1953125, 1,  -0.625, 0.1953125)
loadings <- matrix(loadings, 2, 3, byrow = TRUE)
omxCheckError(mxMatrix("Diag", free=FALSE, values=loadings, name="L", byrow=TRUE),
		paste("Non-square matrix attempted for diagonal matrix constructor",
			  "mxMatrix(\"Diag\", free = FALSE, values = loadings,", 
			  "name = \"L\", byrow = TRUE)"))	
#Error check squareness of a Standardized matrix
loadings <- c(1, -0.625, 0.1953125, 1,  -0.625, 0.1953125)
loadings <- matrix(loadings, 2, 3, byrow = TRUE)
omxCheckError(mxMatrix("Stand", free=FALSE, values=loadings, name="L", byrow=TRUE),
		paste("Non-square matrix attempted for standardized matrix constructor",
			  "mxMatrix(\"Stand\", free = FALSE, values = loadings,", 
			  "name = \"L\", byrow = TRUE)"))		
#Error check squareness of a subdiagonal matrix
loadings <- c(1, -0.625, 0.1953125, 1,  -0.625, 0.1953125)
loadings <- matrix(loadings, 2, 3, byrow = TRUE)
omxCheckError(mxMatrix("Sdiag", free=FALSE, values=loadings, name="L", byrow=TRUE),
		paste("Non-square matrix attempted for subdiagonal matrix constructor",
		      "mxMatrix(\"Sdiag\", free = FALSE, values = loadings,", 
			  "name = \"L\", byrow = TRUE)"))
#Warning Check ignore values warning for Identity Matrix
loadings <- c(1, -0.625, 0.1953125, 1,  -0.625, 0.1953125, 1,  -0.625, 0.1953125)
loadings <- matrix(loadings, 3, 3, byrow = TRUE)
omxCheckWarning(mxMatrix("Iden", free=FALSE, values=loadings, name="L", byrow=TRUE),
        paste("Ignoring values matrix for identity matrix construction",
              "mxMatrix(\"Iden\", free = FALSE, values = loadings,",
              "name = \"L\", byrow = TRUE)"))
#Warning Check ignore values warning for Identity Matrix
loadings <- c(1, -0.625, 0.1953125, 1,  -0.625, 0.1953125, 1,  -0.625, 0.1953125)
loadings <- matrix(loadings, 3, 3, byrow = TRUE)
omxCheckWarning(mxMatrix("Zero", free=FALSE, values=loadings, name="L", byrow=TRUE),
	    paste("Ignoring values matrix for zero matrix constructor",
			  "mxMatrix(\"Zero\", free = FALSE, values = loadings,",
			  "name = \"L\", byrow = TRUE)"))
#Warning Check ignore free matrix for Zero Matrix
omxCheckWarning(mxMatrix("Zero", free = TRUE, name = "L", nrow = 3, ncol = 3, byrow = TRUE),
	    paste("Ignoring free matrix for zero matrix constructor",
			  "mxMatrix(\"Zero\", free = TRUE,",
			  "name = \"L\", nrow = 3, ncol = 3, byrow = TRUE)"))
#Warning Check ignore free matrix for Identity Matrix
omxCheckWarning(mxMatrix("Iden", free = TRUE, name = "L", nrow = 3, ncol = 3, byrow = TRUE),
		paste("Ignoring free matrix for identity matrix construction",
			  "mxMatrix(\"Iden\", free = TRUE,",
			  "name = \"L\", nrow = 3, ncol = 3, byrow = TRUE)"))			
