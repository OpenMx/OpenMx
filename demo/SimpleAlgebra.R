library(OpenMx)

O <- mxMatrix("Full", c(3), nrow=1, ncol=1, name="O")
A <- mxMatrix("Full", c(1, 2), nrow=1, ncol=2, name="A")
B <- mxMatrix("Full", c(3, 4), ncol=1, nrow=2, name="B")
AlgA <- mxAlgebra(O, "AlgA")
AlgAB <- mxAlgebra(A %*% B, "AlgAB")
AlgdoubleMult <- mxAlgebra(A %*% B %*% O, "AlgdoubleMult")

model <- mxModel()
model <- mxModel(model, A)
model <- mxModel(model, O)

# Test 1: Algebra is just a matrix.

model <- mxModel(model, AlgA)
model <- mxModel(model, mxAlgebraObjective("AlgA"))
model <- mxRun(model)

outputA <- model[["AlgA"]]@result
valA <- O[1,1]
diffA <- (valA - outputA[1,1]) / valA
diffA

# Test 2: Algebra is a multiply.

model <- mxModel(model, B)
model <- mxModel(model, AlgAB)
model <- mxModel(model, mxAlgebraObjective("AlgAB"))
model <- mxRun(model)

outputAB <- model[["AlgAB"]]@result
AB <- (A@values %*% B@values)[1,1]
diffAB <- (AB - outputAB[1,1]) / AB
diffAB

# Test 2: Algebra is two multiplies

model <- mxModel(model, AlgdoubleMult)
model <- mxModel(model, mxAlgebraObjective("AlgdoubleMult"))
model <- mxRun(model)

outputdoubleMult <- model[["AlgdoubleMult"]]@result
doubleMult <- ((A@values%*%B@values) %*% O@values)[1,1]
diffdoubleMult <- (doubleMult - outputdoubleMult[1,1]) / doubleMult
diffdoubleMult

c(valA, outputA[1,1], diffA)
c(AB, outputAB[1,1],diffAB)
c(doubleMult, outputdoubleMult[1,1],diffdoubleMult)

