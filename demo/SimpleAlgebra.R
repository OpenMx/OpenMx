library(OpenMx)

modelFit <-function(model, objective) {
	matList <- omxGenerateMatrixList(model)
	algList <- omxGenerateAlgebraList(model)
	varList <- omxGenerateParameterList(model)
	startVals <- c() # omxGenerateValueList(model)
	print(algList)
	return(testFit(objective, startVals, c(), matList, varList, algList))
}

testFit <- function(objective, startVals=c(), bounds=c(), matList=list(), varList=list(), algList=list(), data=c(), state=c()) {
	return(.Call("callNPSOL", objective, startVals, bounds, matList, varList, algList, data, state));
}

# Test 1: Algebra is just a matrix.

O <- mxMatrix("Full", c(1), nrow=1, ncol=1, name="O")
A <- mxMatrix("Full", c(1, 2), nrow=1, ncol=2, name="A")
B <- mxMatrix("Full", c(3, 4), ncol=1, nrow=2, name="B")
AlgA <- mxAlgebra(O, name="AlgA")
AlgAB <- mxAlgebra(A %*% B, name="AlgAB")
AlgRecurse <- mxAlgebra(A %*% B %*% O, name="AlgRecurse")

model <- mxModel()
model <- mxModel(model, A)
model <- mxModel(model, O)

model <- mxModel(model, AlgA)
objectiveA <- mxAlgebraObjective(model, "AlgA")
outputA <- modelFit(model, objectiveA)
outputA

valA <- O[1,1]
diffA <- (valA - outputA$minimum) / valA
diffA

model <- mxModel(model, B)
model <- mxModel(model, AlgAB)
objectiveAB <- mxAlgebraObjective(model, "AlgAB")
outputAB <- modelFit(model, objectiveAB)
outputAB
AB <- (A@values %*% B@values)[1,1]
diffAB <- (AB - outputAB$minimum) / AB
diffAB

model <- mxModel(model, AlgRecurse)
objectiveRecurse <- mxAlgebraObjective(model, "AlgRecurse")
outputRecurse <- modelFit(model, objectiveRecurse)
outputRecurse
recurse <- ((A@values%*%B@values) %*% O@values)[1,1]
diffRecurse <- (recurse - outputRecurse$minimum) / recurse
diffRecurse

c(valA, outputA$minimum, diffA)
c(AB, outputAB$minimum,diffAB)
c(recurse, outputRecurse$minimum,diffRecurse)

