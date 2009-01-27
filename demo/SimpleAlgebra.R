library(OpenMx)

modelFit <-function(model, objective) {
	matList <- omxGenerateMatrixList(model)
	algList <- omxGenerateAlgebraList(model)
	varList <- omxGenerateParameterList(model)
	startVals <- c() # omxGenerateValueList(model)
	print(algList)
	print(objective)
	return(testFit(objective, startVals, c(), matList, varList, algList))
}

testFit <- function(objective, startVals=c(), bounds=c(), matList=list(), varList=list(), algList=list(), data=c(), state=c()) {
	return(.Call("callNPSOL", objective, startVals, bounds, matList, varList, algList, data, state));
}

O <- mxMatrix("Full", c(3), nrow=1, ncol=1, name="O")
A <- mxMatrix("Full", c(1, 2), nrow=1, ncol=2, name="A")
B <- mxMatrix("Full", c(3, 4), ncol=1, nrow=2, name="B")
AlgA <- mxAlgebra(O, name="AlgA")
AlgAB <- mxAlgebra(A %*% B, name="AlgAB")
AlgdoubleMult <- mxAlgebra(A %*% B %*% O, name="AlgdoubleMult")

model <- mxModel()
model <- mxModel(model, A)
model <- mxModel(model, O)

# Test 1: Algebra is just a matrix.

model <- mxModel(model, AlgA)
objectiveA <- mxAlgebraObjective(model, "AlgA")
outputA <- modelFit(model, objectiveA)
outputA

valA <- O[1,1]
diffA <- (valA - outputA$minimum) / valA
diffA

# Test 2: Algebra is a multiply.

model <- mxModel(model, B)
model <- mxModel(model, AlgAB)
objectiveAB <- mxAlgebraObjective(model, "AlgAB")
outputAB <- modelFit(model, objectiveAB)
outputAB
AB <- (A@values %*% B@values)[1,1]
diffAB <- (AB - outputAB$minimum) / AB
diffAB

# Test 2: Algebra is two multiplies

model <- mxModel(model, AlgdoubleMult)
objectivedoubleMult <- mxAlgebraObjective(model, "AlgdoubleMult")
outputdoubleMult <- modelFit(model, objectivedoubleMult)
outputdoubleMult
doubleMult <- ((A@values%*%B@values) %*% O@values)[1,1]
diffdoubleMult <- (doubleMult - outputdoubleMult$minimum) / doubleMult
diffdoubleMult

c(valA, outputA$minimum, diffA)
c(AB, outputAB$minimum,diffAB)
c(doubleMult, outputdoubleMult$minimum,diffdoubleMult)

