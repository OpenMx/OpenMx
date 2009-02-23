library(OpenMx)

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
model <- mxModel(model, mxAlgebraObjective(algebra = "AlgA"))
model <- mxJobRun(model)

outputA <- model@output
valA <- O[1,1]
diffA <- (valA - outputA$minimum) / valA
diffA

# Test 2: Algebra is a multiply.

model <- mxModel(model, B)
model <- mxModel(model, AlgAB)
model <- mxModel(model, mxAlgebraObjective(algebra = "AlgAB"))
model <- mxJobRun(model)

outputAB <- model@output
AB <- (A@values %*% B@values)[1,1]
diffAB <- (AB - outputAB$minimum) / AB
diffAB

# Test 2: Algebra is two multiplies

model <- mxModel(model, AlgdoubleMult)
model <- mxModel(model, mxAlgebraObjective(algebra = "AlgdoubleMult"))
model <- mxJobRun(model)

outputdoubleMult <- model@output
doubleMult <- ((A@values%*%B@values) %*% O@values)[1,1]
diffdoubleMult <- (doubleMult - outputdoubleMult$minimum) / doubleMult
diffdoubleMult

c(valA, outputA$minimum, diffA)
c(AB, outputAB$minimum,diffAB)
c(doubleMult, outputdoubleMult$minimum,diffdoubleMult)

