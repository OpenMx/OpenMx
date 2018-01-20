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

# Error Detection at Runtime: ###

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
	      "A cycle has been detected in model 'model' . It involved the following elements: 'model2.B' and 'model1.A'
A common trigger for this error is not providing a name string as the first parameter to mxModel.")

A <- mxMatrix('Full', 1, 1, labels = 'B[1,1]', name = 'A')
B <- mxMatrix('Full', 1, 1, labels = 'C[1,1]', name = 'B')
C <- mxMatrix('Full', 1, 1, labels = 'A[1,1]', name = 'C')
model <- mxModel('model', A, B, C)
omxCheckError(mxRun(model),
	      "A cycle has been detected in model 'model' . It involved the following elements: 'B', 'C', and 'A'
A common trigger for this error is not providing a name string as the first parameter to mxModel.")

A <- mxMatrix('Full', 2, 2, name = 'A')
B <- mxMatrix('Full', 2, 2, name = 'B')
C <- mxAlgebra(A + B, name = 'C')
D <- mxMatrix('Full', 1, 1, labels = 'C[2,3]', name = 'D')
model <- mxModel('model', A, B, C, D,
                 mxFitFunctionAlgebra("D"))
omxCheckError(mxRun(model), 
	"Requested improper value (2, 3) from (2, 2) matrix 'model.C'")

# Error Detection at Construction, by argument to function mxMatrix():  ###

#Error check for type:
loadings <- c(1, -0.625, 0.1953125, 1,  -0.625, 0.1953125)
loadings <- matrix(loadings, 2, 3, byrow = TRUE)
omxCheckError(
  mxMatrix("wrong", free=FALSE, values=loadings, name="L", byrow=TRUE),
  '\'type\' must be one of: Diag Full Iden Lower Stand Sdiag Symm Unit Zero in mxMatrix("wrong", free = FALSE, values = loadings, name = "L", byrow = TRUE)')
#Error check when both nrow and ncol are NA:
loadings <- c(1, -0.625, 0.1953125, 1,  -0.625, 0.1953125)
omxCheckError(mxMatrix("Full", free=FALSE, values=loadings, name="L", byrow=TRUE),
              paste("you must specify 'nrow' and 'ncol' arguments in mxMatrix(\"Full\",",
                    "free = FALSE, values = loadings, name = \"L\", byrow = TRUE)"))
#Error check when nrow or ncol are non-scalar or non-numeric:
omxCheckError(
  mxMatrix("Full", free=FALSE, values=loadings, name="L", nrow=c(2,9), ncol=3, byrow=TRUE),
  '\'nrow\' argument to mxMatrix function must be either NA or a single numeric value in mxMatrix("Full", free = FALSE, values = loadings, name = "L", nrow = c(2, 9), ncol = 3, byrow = TRUE)')
omxCheckError(
  mxMatrix("Full", free=FALSE, values=loadings, name="L", nrow=2, ncol="3", byrow=TRUE),
  '\'ncol\' argument to mxMatrix function must be either NA or a single numeric value in mxMatrix("Full", free = FALSE, values = loadings, name = "L", nrow = 2, ncol = "3", byrow = TRUE)')
#Error check for 'free' argument:
loadings <- c(1, -0.625, 0.1953125, 1,  -0.625, 0.1953125)
loadings <- matrix(loadings, 2, 3, byrow = TRUE)
omxCheckError(mxMatrix("Full", free="FALSE", values=loadings, name="L", byrow=TRUE),
              paste("'free' argument to mxMatrix function must be of logical type",
                    "in mxMatrix(\"Full\", free = \"FALSE\", values = loadings,",
                    "name = \"L\", byrow = TRUE)"))
omxCheckError(
  mxMatrix("Full", free=NA, values=loadings, name="L", byrow=TRUE),
  '\'free\' argument to mxMatrix function cannot contain an NA in mxMatrix("Full", free = NA, values = loadings, name = "L", byrow = TRUE)')
#Error check for values: 
loadings <- c(1, -0.625, 0.1953125, 1,  "h", 0.1953125)
loadings <- matrix(loadings, 2, 3, byrow = TRUE)
omxCheckError(mxMatrix("Full", free=FALSE, values=loadings, name="L", byrow=TRUE),
        paste("'values' argument to mxMatrix function must be of numeric type",
              "in mxMatrix(\"Full\", free = FALSE, values = loadings,",
              "name = \"L\", byrow = TRUE)"))
#Error check for labels:
loadings <- c(1, -0.625, 0.1953125, 1,  -0.625, 0.1953125)
omxCheckError(mxMatrix("Full", free=TRUE, values=matrix(0,2,3), labels=loadings, name="L", byrow=TRUE),
              "'labels' argument to mxMatrix function must be of character type in mxMatrix(\"Full\", free = TRUE, values = matrix(0, 2, 3), labels = loadings, name = \"L\", byrow = TRUE)")
#Error check for lbound and ubound:
loadings <- c(1, -0.625, 0.1953125, 1,  "h", 0.1953125)
loadings <- matrix(loadings, 2, 3, byrow = TRUE)
omxCheckError(
  mxMatrix("Full", free=TRUE, values=matrix(0,2,3), lbound=loadings, name="L", byrow=TRUE),
  "'lbound' argument to mxMatrix function must be of numeric type in mxMatrix(\"Full\", free = TRUE, values = matrix(0, 2, 3), lbound = loadings, name = \"L\", byrow = TRUE)")
omxCheckError(
  mxMatrix("Full", free=TRUE, values=matrix(0,2,3), ubound=loadings, name="L", byrow=TRUE),
  "'ubound' argument to mxMatrix function must be of numeric type in mxMatrix(\"Full\", free = TRUE, values = matrix(0, 2, 3), ubound = loadings, name = \"L\", byrow = TRUE)")
#Error check for name:
omxCheckError(
  mxMatrix("Full", free=FALSE, values=loadings, nrow=2, ncol=3, name=TRUE, byrow=TRUE),
  "'name' argument must be a character vector in mxMatrix(\"Full\", free = FALSE, values = loadings, nrow = 2, ncol = 3, name = TRUE, byrow = TRUE)"
)
#Error check for condenseSlots:
omxCheckError(
  mxMatrix("Full", free=FALSE, values=matrix(0,2,3), name="L", byrow=TRUE, condenseSlots=NA),
  '\'condenseSlots\' argument to mxMatrix function cannot contain an NA in mxMatrix("Full", free = FALSE, values = matrix(0, 2, 3), name = "L", byrow = TRUE, condenseSlots = NA)')
#Error check for S4 class supplied as argument to mxMatrix():
omxCheckError(
  mxMatrix("Full", free=FALSE, values=mxMatrix("Iden",nrow=2,name="I")),
  '\'values\' argument to mxMatrix function must be a scalar, a vector, or a matrix in mxMatrix("Full", free = FALSE, values = mxMatrix("Iden", nrow = 2, name = "I"))')
#Error check for conflicting dimensions of inputs:
omxCheckError(
  mxMatrix("Full",free=matrix(F,2,3),values=matrix(0,3,2),name="L"),
  'Two or more matrix inputs have different dimensions. Use the \'nrow\' and \'ncol\' arguments in mxMatrix("Full", free = matrix(F, 2, 3), values = matrix(0, 3, 2), name = "L")')
omxCheckError(
  mxMatrix("Full",free=matrix(F,2,3),values=matrix(0,3,2),nrow=2,name="L"),
  'Two or more matrix inputs have different dimensions. Use the \'nrow\' and \'ncol\' arguments in mxMatrix("Full", free = matrix(F, 2, 3), values = matrix(0, 3, 2), nrow = 2, name = "L")')
              


# Type-specific errors and warnings: ###

#Error check squareness of a Diagonal MxMatrix:
loadings <- c(1, -0.625, 0.1953125, 1,  -0.625, 0.1953125)
loadings <- matrix(loadings, 2, 3, byrow = TRUE)
omxCheckError(mxMatrix("Diag", free=FALSE, values=loadings, name="L", byrow=TRUE),
              paste("non-square MxMatrix attempted in 'nrow' and 'ncol' arguments to",
                    "mxMatrix(\"Diag\", free = FALSE, values = loadings,", 
                    "name = \"L\", byrow = TRUE)"))	
#Error check diagonality of values argument:
omxCheckError(
  mxMatrix("Diag", free=FALSE, values=matrix(1,3,3), name="L", byrow=TRUE),
  '\'values\' matrix of L is not a diagonal matrix in mxMatrix("Diag", free = FALSE, values = matrix(1, 3, 3), name = "L", byrow = TRUE)')
#Error check diagonality of 'free' argument:
omxCheckError(
  mxMatrix("Diag", free=matrix(TRUE,3,3), values=diag(3,3), name="L", byrow=TRUE),
  '\'free\' matrix of L has TRUE on non-diagonal in mxMatrix("Diag", free = matrix(TRUE, 3, 3), values = diag(3, 3), name = "L", byrow = TRUE)')

#Error check for attempted non-square dimensions for Identity MxMatrix:
omxCheckError(
  mxMatrix("Iden",nrow=3,ncol=4,name="I"),
  'non-square MxMatrix attempted for Identity MxMatrix construction in mxMatrix("Iden", nrow = 3, ncol = 4, name = "I")')
#Warning Check ignore values warning for Identity MxMatrix:
loadings <- c(1, -0.625, 0.1953125, 1,  -0.625, 0.1953125, 1,  -0.625, 0.1953125)
loadings <- matrix(loadings, 3, 3, byrow = TRUE)
omxCheckWarning(mxMatrix("Iden", free=FALSE, values=loadings, name="I", byrow=TRUE),
                paste("ignoring 'values' matrix for Identity MxMatrix construction in",
                      "mxMatrix(\"Iden\", free = FALSE, values = loadings,",
                      "name = \"I\", byrow = TRUE)"))
#Warning Check ignore free matrix for Identity MxMatrix:
omxCheckWarning(mxMatrix("Iden", free = TRUE, name = "I", nrow = 3, ncol = 3, byrow = TRUE),
                paste("ignoring 'free' matrix for Identity MxMatrix construction in",
                      "mxMatrix(\"Iden\", free = TRUE,",
                      "name = \"I\", nrow = 3, ncol = 3, byrow = TRUE)"))
#Warning Check ignore labels matrix for Identity MxMatrix:
omxCheckWarning(
  mxMatrix("Iden", labels=matrix(letters[1:9],3,3), name = "I", nrow = 3, ncol = 3, byrow = TRUE),
  'ignoring \'labels\' matrix for Identity MxMatrix construction in mxMatrix("Iden", labels = matrix(letters[1:9], 3, 3), name = "I", nrow = 3, ncol = 3, byrow = TRUE)')
#Warning Check ignore lbound matrix for Identity MxMatrix:
omxCheckWarning(
  mxMatrix("Iden", lbound=-1, name = "I", nrow = 3, ncol = 3, byrow = TRUE),
  'ignoring \'lbound\' matrix for Identity MxMatrix construction in mxMatrix("Iden", lbound = -1, name = "I", nrow = 3, ncol = 3, byrow = TRUE)')
#Warning Check ignore ubound matrix for Identity MxMatrix:
omxCheckWarning(
  mxMatrix("Iden", ubound=2, name = "I", nrow = 3, ncol = 3, byrow = TRUE),
  'ignoring \'ubound\' matrix for Identity MxMatrix construction in mxMatrix("Iden", ubound = 2, name = "I", nrow = 3, ncol = 3, byrow = TRUE)')


#Error check upper values of a Lower MxMatrix:
loadings <- c(1, -0.625, 0.1953125, 1,  -0.625, 0.1953125, 1,  -0.625, 0.1953125)
loadings <- matrix(loadings, 3, 3, byrow = TRUE)
omxCheckError(mxMatrix("Lower", free=FALSE, values=loadings, name="L", byrow=TRUE),
        paste("upper triangle of 'values' matrix in Lower MxMatrix 'L' is not all zeros in",
              "mxMatrix(\"Lower\", free = FALSE, values = loadings,", 
              "name = \"L\", byrow = TRUE)"))
#Error check free status of upper values of a Lower MxMatrix:
omxCheckError(
  mxMatrix("Lower",free=matrix(T,3,3),nrow=3,ncol=3,name="L",byrow=T),
  'upper triangle of \'free\' matrix in Lower MxMatrix \'L\' is not all fixed in mxMatrix("Lower", free = matrix(T, 3, 3), nrow = 3, ncol = 3, name = "L", byrow = T)')
#Error check labels of upper triangle of a Lower MxMatrix:
omxCheckError(
  mxMatrix("Lower",labels=matrix(letters[1:9],3,3),nrow=3,ncol=3,name="L",byrow=T),
  'upper triangle of \'labels\' matrix in Lower MxMatrix \'L\' is not all NAs in mxMatrix("Lower", labels = matrix(letters[1:9], 3, 3), nrow = 3, ncol = 3, name = "L", byrow = T)')
#Error check lbound of upper triangle of a Lower MxMatrix:
omxCheckError(
  mxMatrix("Lower",lbound=matrix(-1,3,3),nrow=3,ncol=3,name="L",byrow=T),
  'upper triangle of \'lbound\' matrix in Lower MxMatrix \'L\' is not all NAs in mxMatrix("Lower", lbound = matrix(-1, 3, 3), nrow = 3, ncol = 3, name = "L", byrow = T)')
#Error check ubound of upper triangle of a Lower MxMatrix:
omxCheckError(
  mxMatrix("Lower",ubound=matrix(1,3,3),nrow=3,ncol=3,name="L",byrow=T),
  'upper triangle of \'ubound\' matrix in Lower MxMatrix \'L\' is not all NAs in mxMatrix("Lower", ubound = matrix(1, 3, 3), nrow = 3, ncol = 3, name = "L", byrow = T)')
#Error check squareness of a Lower MxMatrix:
loadings <- c(1, -0.625, 0.1953125, 1,  -0.625, 0.1953125)
loadings <- matrix(loadings, 2, 3, byrow = TRUE)
omxCheckError(mxMatrix("Lower", free=FALSE, values=loadings, name="L", byrow=TRUE),
		paste("non-square MxMatrix attempted in 'nrow' and 'ncol' arguments to",
			  "mxMatrix(\"Lower\", free = FALSE, values = loadings,", 
			  "name = \"L\", byrow = TRUE)"))
#Error check bad number of elements supplied to Lower MxMatrix:
omxCheckError(
  mxMatrix("Lower",nrow=3,values=1:7,name="L"),
  'illegal number of elements (7) for \'values\' matrix in Lower MxMatrix construction mxMatrix("Lower", nrow = 3, values = 1:7, name = "L")')

#Error check squareness of a Subdiagonal MxMatrix:
loadings <- c(1, -0.625, 0.1953125, 1,  -0.625, 0.1953125)
loadings <- matrix(loadings, 2, 3, byrow = TRUE)
omxCheckError(mxMatrix("Sdiag", free=FALSE, values=loadings, name="L", byrow=TRUE),
              paste("non-square MxMatrix attempted in 'nrow' and 'ncol' arguments to",
                    "mxMatrix(\"Sdiag\", free = FALSE, values = loadings,", 
                    "name = \"L\", byrow = TRUE)"))
#Error check bad number of elements supplied to Subdiagonal MxMatrix:
omxCheckError(
  mxMatrix("Sdiag",nrow=3,ncol=3,values=1:4,name="S"),
'illegal number of elements (4) for \'values\' matrix of Subdiagonal MxMatrix constructionmxMatrix("Sdiag", nrow = 3, ncol = 3, values = 1:4, name = "S")')
#Error check free status of upper values of a Subdiagonal MxMatrix:
omxCheckError(
  mxMatrix("Sdiag",free=matrix(T,3,3),nrow=3,ncol=3,name="S",byrow=T),
  'upper triangle or diagonal of \'free\' matrix in Subdiagonal MxMatrix \'S\' is not all fixed in mxMatrix("Sdiag", free = matrix(T, 3, 3), nrow = 3, ncol = 3, name = "S", byrow = T)')
#Error check labels of Subdiagonal MxMatrix:
omxCheckError(
  mxMatrix("Sdiag",labels=matrix(letters[1:9],3,3),nrow=3,ncol=3,name="S",byrow=T),
  'upper triangle or diagonal of \'labels\' matrix in Subdiagonal MxMatrix \'S\' is not all NAs in mxMatrix("Sdiag", labels = matrix(letters[1:9], 3, 3), nrow = 3, ncol = 3, name = "S", byrow = T)')
#Error check lbound of Subdiagonal MxMatrix:
omxCheckError(
  mxMatrix("Sdiag",lbound=matrix(-1,3,3),nrow=3,ncol=3,name="S",byrow=T),
  'upper triangle or diagonal of \'lbound\' matrix in Subdiagonal MxMatrix \'S\' is not all NAs in mxMatrix("Sdiag", lbound = matrix(-1, 3, 3), nrow = 3, ncol = 3, name = "S", byrow = T)')
#Error check ubound Subdiagonal MxMatrix:
omxCheckError(
  mxMatrix("Sdiag",ubound=matrix(1,3,3),nrow=3,ncol=3,name="S",byrow=T),
  'upper triangle or diagonal of \'ubound\' matrix in Subdiagonal MxMatrix \'S\' is not all NAs in mxMatrix("Sdiag", ubound = matrix(1, 3, 3), nrow = 3, ncol = 3, name = "S", byrow = T)')

#Error check squareness of a Standardized MxMatrix:
loadings <- c(1, -0.625, 0.1953125, 1,  -0.625, 0.1953125)
loadings <- matrix(loadings, 2, 3, byrow = TRUE)
omxCheckError(mxMatrix("Stand", free=FALSE, values=loadings, name="S", byrow=TRUE),
              paste("non-square MxMatrix attempted in 'nrow' and 'ncol' arguments to",
                    "mxMatrix(\"Stand\", free = FALSE, values = loadings,", 
                    "name = \"S\", byrow = TRUE)"))
#Error check bad number of elements supplied to Standardized MxMatrix:
omxCheckError(
  mxMatrix("Stand",nrow=3,ncol=3,values=1:4,name="S"),
  'illegal number of elements (4) for \'values\' matrix in Standardized MxMatrix construction mxMatrix("Stand", nrow = 3, ncol = 3, values = 1:4, name = "S")')
#Error check symmetry of Standardized MxMatrix:
omxCheckError(
  mxMatrix("Stand",values=matrix(c(1,2,3,1),2,2),name="S"),
  '\'values\' matrix of Standardized MxMatrix \'S\' is not symmetric in mxMatrix("Stand", values = matrix(c(1, 2, 3, 1), 2, 2), name = "S")')
#Error check diagonal of 1s in Standardized MxMatrix:
omxCheckError(
  mxMatrix("Stand",values=matrix(c(1,2,2,2),2,2),name="S"),
  "'values' matrix of Standardized MxMatrix 'S' is not 1s along the diagonal in mxMatrix(\"Stand\", values = matrix(c(1, 2, 2, 2), 2, 2), name = \"S\")")
omxCheckError(
  mxMatrix("Stand",values=matrix(c(1-1e-15,2,2,1-1e-15),2,2),name="S"),
  "'values' matrix of Standardized MxMatrix 'S' is very near, but not equal to, 1s along the diagonal in mxMatrix(\"Stand\", values = matrix(c(1 - 1e-15, 2, 2, 1 - 1e-15), 2, 2), name = \"S\")")
#Error check 'free' argument to Standardized MxMatrix:
omxCheckError(
  mxMatrix("Stand",nrow=2,free=matrix(c(F,F,T,F),2,2),name="S"),
  "'free' matrix of Standardized MxMatrix 'S' is not symmetric in mxMatrix(\"Stand\", nrow = 2, free = matrix(c(F, F, T, F), 2, 2), name = \"S\")")
omxCheckError(
  mxMatrix("Stand",nrow=2,free=matrix(c(T,F,F,F),2,2),name="S"),
  "'free' matrix of Standardized MxMatrix 'S' is not fixed along the diagonal in mxMatrix(\"Stand\", nrow = 2, free = matrix(c(T, F, F, F), 2, 2), name = \"S\")")
#Error check 'labels' argument to Standardized MxMatrix:
omxCheckError(
  mxMatrix("Stand",nrow=2,labels=c("a","b","b","c"),name="S"),
  '\'labels\' matrix of Standardized MxMatrix \'S\' is not NA along the diagonal in mxMatrix("Stand", nrow = 2, labels = c("a", "b", "b", "c"), name = "S")')
omxCheckError(
  mxMatrix("Stand",nrow=2,labels=c("a","b","c","d"),name="S"),
  '\'labels\' matrix of Standardized MxMatrix \'S\' is not symmetric in mxMatrix("Stand", nrow = 2, labels = c("a", "b", "c", "d"), name = "S")')
#Error check 'lbound' argument to Standardized MxMatrix:
omxCheckError(
  mxMatrix("Stand",nrow=2,lbound=c(NA,-1,-2,NA),name="S"),
  "'lbound' matrix of Standardized MxMatrix 'S' is not symmetric in mxMatrix(\"Stand\", nrow = 2, lbound = c(NA, -1, -2, NA), name = \"S\")")
omxCheckError(
  mxMatrix("Stand",nrow=2,lbound=rep(-1,4),name="S"),
  "'lbound' matrix of Standardized MxMatrix 'S' is not NA along the diagonal in mxMatrix(\"Stand\", nrow = 2, lbound = rep(-1, 4), name = \"S\")")
#Error check 'ubound' argument to Standardized MxMatrix:
omxCheckError(
  mxMatrix("Stand",nrow=2,ubound=c(NA,-1,-2,NA),name="S"),
  "'ubound' matrix of Standardized MxMatrix 'S' is not symmetric in mxMatrix(\"Stand\", nrow = 2, ubound = c(NA, -1, -2, NA), name = \"S\")")
omxCheckError(
  mxMatrix("Stand",nrow=2,ubound=rep(-1,4),name="S"),
  "'ubound' matrix of Standardized MxMatrix 'S' is not NA along the diagonal in mxMatrix(\"Stand\", nrow = 2, ubound = rep(-1, 4), name = \"S\")")

#Error check squareness of a Symmetrical MxMatrix:
loadings <- c(1, -0.625, 0.1953125, 1,  -0.625, 0.1953125)
loadings <- matrix(loadings, 2, 3, byrow = TRUE)
omxCheckError(mxMatrix("Symm", free=FALSE, values=loadings, name="L", byrow=TRUE),
	    paste("non-square MxMatrix attempted in 'nrow' and 'ncol' arguments to",
			  "mxMatrix(\"Symm\", free = FALSE, values = loadings,", 
			  "name = \"L\", byrow = TRUE)"))
#Error check bad number of elements supplied to Symmetrical MxMatrix:
omxCheckError(
  mxMatrix("Symm",nrow=3,values=1:10,name="S"),
  "illegal number of elements (10) for 'values' matrix in Symmetric MxMatrix constructionmxMatrix(\"Symm\", nrow = 3, values = 1:10, name = \"S\")")
#Error check 'values' argument to Symmetrical MxMatrix:
omxCheckError(
  mxMatrix("Symm",values=c(3,4,NA,5),nrow=2,name="S"),
  "NAs in 'values' matrix of Symmetric MxMatrix 'S' are not symmetric in mxMatrix(\"Symm\", values = c(3, 4, NA, 5), nrow = 2, name = \"S\")")
omxCheckError(
  mxMatrix("Symm",values=c(1,2,3,4),nrow=2,name="S"),
  "'values' matrix of Symmetric MxMatrix 'S' is not symmetric in mxMatrix(\"Symm\", values = c(1, 2, 3, 4), nrow = 2, name = \"S\")")
#Error check 'labels' argument to Symmetrical MxMatrix:
omxCheckError(
  mxMatrix("Symm",values=diag(2),labels=c("a","b","c","d"),nrow=2,name="S"),
  '\'labels\' matrix of Symmetric MxMatrix \'S\' is not symmetric in mxMatrix("Symm", values = diag(2), labels = c("a", "b", "c", "d"), nrow = 2, name = "S")')
#Error check 'lbound' argument to Symmetrical MxMatrix:
omxCheckError(
  mxMatrix("Symm",nrow=2,lbound=c(1,1,2,1),name="S"),
  "'lbound' matrix of Symmetric MxMatrix 'S' is not symmetric in mxMatrix(\"Symm\", nrow = 2, lbound = c(1, 1, 2, 1), name = \"S\")")
#Error check 'ubound' argument to Symmetrical MxMatrix:
omxCheckError(
  mxMatrix("Symm",nrow=2,ubound=c(1,1,2,1),name="S"),
  "'ubound' matrix of Symmetric MxMatrix 'S' is not symmetric in mxMatrix(\"Symm\", nrow = 2, ubound = c(1, 1, 2, 1), name = \"S\")")

#Warning check 'values' argument to Unit MxMatrix:
omxCheckWarning(
  mxMatrix("Unit",nrow=2,ncol=3,values=1:6,name="U"),
  'ignoring \'values\' matrix for Unit MxMatrix construction in mxMatrix("Unit", nrow = 2, ncol = 3, values = 1:6, name = "U")')
#Warning check 'labels' argument to Unit MxMatrix:
omxCheckWarning(
  mxMatrix("Unit",nrow=2,ncol=3,labels=letters[1:6],name="U"),
  'ignoring \'labels\' matrix for Unit MxMatrix construction in mxMatrix("Unit", nrow = 2, ncol = 3, labels = letters[1:6], name = "U")')
#Warning check 'free' argument to Unit MxMatrix:
omxCheckWarning(
  mxMatrix("Unit",nrow=2,ncol=3,free=T,name="U"),
  'ignoring \'free\' matrix for Unit MxMatrix construction in mxMatrix("Unit", nrow = 2, ncol = 3, free = T, name = "U")')
#Warning check 'lbound' argument to Unit MxMatrix:
omxCheckWarning(
  mxMatrix("Unit",nrow=2,ncol=3,lbound=0,name="U"),
  'ignoring \'lbound\' matrix for Unit MxMatrix construction in mxMatrix("Unit", nrow = 2, ncol = 3, lbound = 0, name = "U")')
#Warning check 'ubound' argument to Unit MxMatrix:
omxCheckWarning(
  mxMatrix("Unit",nrow=2,ncol=3,ubound=2,name="U"),
  'ignoring \'ubound\' matrix for Unit MxMatrix construction in mxMatrix("Unit", nrow = 2, ncol = 3, ubound = 2, name = "U")')

#Warning Check ignore values warning for Zero MxMatrix:
loadings <- c(1, -0.625, 0.1953125, 1,  -0.625, 0.1953125, 1,  -0.625, 0.1953125)
loadings <- matrix(loadings, 3, 3, byrow = TRUE)
omxCheckWarning(mxMatrix("Zero", free=FALSE, values=loadings, name="L", byrow=TRUE),
	    paste("ignoring 'values' matrix for Zero MxMatrix construction in",
			  "mxMatrix(\"Zero\", free = FALSE, values = loadings,",
			  "name = \"L\", byrow = TRUE)"))
#Warning Check ignore free matrix for Zero MxMatrix:
omxCheckWarning(mxMatrix("Zero", free = TRUE, name = "L", nrow = 3, ncol = 3, byrow = TRUE),
	    paste("ignoring 'free' matrix for Zero MxMatrix construction in",
			  "mxMatrix(\"Zero\", free = TRUE,",
			  "name = \"L\", nrow = 3, ncol = 3, byrow = TRUE)"))
#Warning check 'labels' argument to Zero MxMatrix:
omxCheckWarning(
  mxMatrix("Zero",nrow=2,ncol=3,labels=letters[1:6],name="Z"),
  'ignoring \'labels\' matrix for Zero MxMatrix construction in mxMatrix("Zero", nrow = 2, ncol = 3, labels = letters[1:6], name = "Z")')
#Warning check 'lbound' argument to Zero MxMatrix:
omxCheckWarning(
  mxMatrix("Zero",nrow=2,ncol=3,lbound=0,name="Z"),
  'ignoring \'lbound\' matrix for Zero MxMatrix construction in mxMatrix("Zero", nrow = 2, ncol = 3, lbound = 0, name = "Z")')
#Warning check 'ubound' argument to Zero MxMatrix:
omxCheckWarning(
  mxMatrix("Zero",nrow=2,ncol=3,ubound=2,name="Z"),
  'ignoring \'ubound\' matrix for Zero MxMatrix construction in mxMatrix("Zero", nrow = 2, ncol = 3, ubound = 2, name = "Z")')


# Test for MxMatrix type retention (when possible) after taking submatrix: ###

D <- mxMatrix("Diag",nrow=3,values=1:3,labels=letters[1:3],free=T,lbound=0,ubound=4,name="D")
omxCheckEquals(class(D[1:2,1:2])[[1]],"DiagMatrix")
omxCheckEquals(class(D[2:3,2:3])[[1]],"DiagMatrix")
omxCheckEquals(class(D[-2,-2])[[1]],"DiagMatrix")
omxCheckEquals(class(D[3,1])[[1]],"DiagMatrix")


I <- mxMatrix("Iden",nrow=3,name="I")
omxCheckEquals(class(I[1:2,1:2])[[1]],"IdenMatrix")
omxCheckEquals(class(I[1,2])[[1]],"IdenMatrix") #<--This passes, which is a bit counterintuitive.

L <- mxMatrix("Lower",nrow=3,values=1:6,name="L")
omxCheckEquals(class(L[1:2,1:2])[[1]],"LowerMatrix")
omxCheckEquals(class(L[1:2,2:3])[[1]],"LowerMatrix")
omxCheckEquals(class(L[2:3,2:3])[[1]],"LowerMatrix")
omxCheckEquals(class(L[1,3])[[1]],"LowerMatrix")
omxCheckEquals(class(L[-2,-2])[[1]],"LowerMatrix")

S <- mxMatrix("Sdiag",nrow=3,values=1:3,name="S")
omxCheckEquals(class(S[1:2,1:2])[[1]],"SdiagMatrix")
omxCheckEquals(class(S[1:2,2:3])[[1]],"SdiagMatrix")
omxCheckEquals(class(S[2:3,2:3])[[1]],"SdiagMatrix")
omxCheckEquals(class(S[-2,-2])[[1]],"SdiagMatrix")

S <- mxMatrix("Stand",nrow=3,values=1:3,name="S")
omxCheckEquals(class(S[1:2,1:2])[[1]],"StandMatrix")
omxCheckEquals(class(S[2:3,2:3])[[1]],"StandMatrix")
omxCheckEquals(class(S[-2,-2])[[1]],"StandMatrix")

S <- mxMatrix("Symm",nrow=3,values=1:6,name="S")
omxCheckEquals(class(S[1:2,1:2])[[1]],"SymmMatrix")
omxCheckEquals(class(S[2:3,2:3])[[1]],"SymmMatrix")
omxCheckEquals(class(S[-2,-2])[[1]],"SymmMatrix")

U <- mxMatrix("Unit",nrow=2,ncol=3,name="U")
omxCheckEquals(class(U[1,2:3])[[1]],"UnitMatrix")

Z <- mxMatrix("Zero",nrow=2,ncol=3,name="Z")
omxCheckEquals(class(Z[2,1:2])[[1]],"ZeroMatrix")


#Modifications to MxMatrices, to be caught as errors at runtime: ###

J <- mxMatrix("Iden",nrow=3,name="I")
J@values[1,2] <- 2
omxCheckError(
  mxRun(mxModel("asdf",J)),
  "'values' matrix of Identity MxMatrix 'I' is not the identity matrix in \"(oops) could not find function mxMatrix\"")

J <- mxMatrix("Iden",nrow=3,name="I",condenseSlots=FALSE)
J@free[1,2] <- TRUE
omxCheckError(
  mxRun(mxModel("asdf",J)),
"'free' matrix of Identity MxMatrix 'I' has a free parameter in \"(oops) could not find function mxMatrix\"")

J <- mxMatrix("Iden",nrow=3,name="I")
J$labels <- matrix(letters[1:9],3,3)
m <- mxModel("asdf",J)
m <- omxSetParameters(model = m, labels = "a", free = T)
omxCheckError(
  mxRun(m),
  "'free' matrix of Identity MxMatrix 'I' has a free parameter in \"(oops) could not find function mxMatrix\"")

U <- mxMatrix("Unit",nrow=2,ncol=3,name="U")
U@values[1,2] <- 2
omxCheckError(
  mxRun(mxModel("asdf",U)),
  "'values' matrix of Unit MxMatrix 'U' has non unit entries in \"(oops) could not find function mxMatrix\"")

U <- mxMatrix("Unit",nrow=2,ncol=3,name="U",condenseSlots=FALSE)
U@free[1,2] <- TRUE
omxCheckError(
  mxRun(mxModel("asdf",U)),
"'free' matrix of Unit MxMatrix 'U' has a free parameter in \"(oops) could not find function mxMatrix\"")

U <- mxMatrix("Unit",nrow=2,ncol=3,name="U")
U$labels <- matrix(letters[1:6],2,3)
m <- mxModel("asdf",U)
m <- omxSetParameters(model = m, labels = "a", free = T)
omxCheckError(
  mxRun(m),
  "'free' matrix of Unit MxMatrix 'U' has a free parameter in \"(oops) could not find function mxMatrix\"")

Z <- mxMatrix("Zero",nrow=2,ncol=3,name="Z")
Z@values[1,2] <- 2
omxCheckError(
  mxRun(mxModel("asdf",Z)),
  "'values' matrix of Zero MxMatrix 'Z' has non-zero entries in \"(oops) could not find function mxMatrix\"")

Z <- mxMatrix("Zero",nrow=2,ncol=3,name="Z",condenseSlots=FALSE)
Z@free[1,2] <- TRUE
omxCheckError(
  mxRun(mxModel("asdf",Z)),
  "'free' matrix of Zero MxMatrix 'Z' has a free parameter in \"(oops) could not find function mxMatrix\"")

Z <- mxMatrix("Zero",nrow=2,ncol=3,name="Z")
Z$labels <- matrix(letters[1:6],2,3)
m <- mxModel("asdf",Z)
m <- omxSetParameters(model = m, labels = "a", free = T)
omxCheckError(
  mxRun(m),
  "'free' matrix of Zero MxMatrix 'Z' has a free parameter in \"(oops) could not find function mxMatrix\"")

Fu <- mxMatrix(type = "Full", nrow = 2, ncol=3, free = rep(F,6), values = 1:6, name="Fu", condenseSlots = T,
               dimnames = list(c("r1","r2"),c("c1","c2","c3")))
Fu@free <- matrix(T,2,3)
omxCheckError(
  mxRun(mxModel("asdf",Fu)),
  "'labels' and 'free' matrices of 'Fu' have different dimensions")

Fu <- mxMatrix(type = "Full", nrow = 2, ncol=3, free = rep(F,6), values = 1:6, name="Fu", condenseSlots = T,
               dimnames = list(c("r1","r2"),c("c1","c2","c3")))
Fu@free <- matrix(T,1,1)
omxCheckError(
  mxRun(mxModel("asdf",Fu)),
  "'values' and 'free' matrices of 'Fu' have different dimensions")

Fu <- mxMatrix(type = "Full", nrow = 2, ncol=3, free = rep(F,6), values = 1:6, name="Fu", condenseSlots = T,
               dimnames = list(c("r1","r2"),c("c1","c2","c3")))
Fu@.condenseSlots <- FALSE
omxCheckError(
  mxRun(mxModel("asdf",Fu)),
  "'labels' and 'values' matrices of 'Fu' have different dimensions")
Fu <- imxConDecMatrixSlots(Fu)
mxRun(mxModel("asdf",Fu)) #<--Should NOT produce an error

Fu <- mxMatrix(type = "Full", nrow = 2, ncol=3, free = rep(F,6), values = 1:6, name="Fu", condenseSlots = T,
               dimnames = list(c("r1","r2"),c("c1","c2","c3")))
Fu@labels <- matrix("nana",1,1)
omxCheckError(
  mxRun(mxModel("asdf",Fu)),
  "'labels' and 'values' matrices of 'Fu' have different dimensions")


#Error and warning checks for replaceMethods: ###

Fu <- mxMatrix(type = "Full", nrow = 2, ncol=3, free = rep(F,6), values = 1:6, name="Fu", condenseSlots = T,
               dimnames = list(c("r1","r2"),c("c1","c2","c3")))
omxCheckWarning(
  Fu$values <- 1:7,
  "data length 7 is not a sub-multiple or multiple of the number of rows 2 for argument 'values' in \"(oops) could not find function mxMatrix\"")
omxCheckWarning(
  Fu$values <- 1:8,
  "data length 8 is not a sub-multiple or multiple of the number of columns 3 for argument 'values' in \"(oops) could not find function mxMatrix\"")
omxCheckError(
  Fu[1:2,1:2] <- "qwer",
  "right-hand side must be MxMatrix object")
omxCheckError(
  dimnames(Fu) <- list(c("a","b")),
  "the 'dimnames' argument to MxMatrix 'Fu' must have a length of 2")
omxCheckError(
  dimnames(Fu) <- c("a","b"),
  "the MxMatrix object 'Fu' has specified dimnames with dimensions 1 x 1 but the matrix is of dimensions 2 x 3")

# test square brackets
m1 <- mxModel("sb", Fu)
v <- 2.5
names(v) <- "sb.Fu[1,1]"
m1 <- omxSetParameters(m1, names(v), values=v, strict=FALSE)
omxCheckEquals(m1$Fu$values[1,1], 2.5)
