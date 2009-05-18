library(OpenMx)

O <- mxMatrix("Full", c(3), nrow=1, ncol=1, name="O")
A <- mxMatrix("Full", c(1, 2), nrow=1, ncol=2, name="A")
B <- mxMatrix("Full", c(3, 4), ncol=1, nrow=2, free = c(TRUE, FALSE), name="B")
algebra <- mxAlgebra(A %*% B - O, "algebra")
algSquared <- mxAlgebra(algebra %*% algebra, "algebraSquared")

model <- mxModel()
model <- mxModel(model, A, B, O)
model <- mxModel(model, algebra, algSquared)
model <- mxModel(model, mxAlgebraObjective("algebraSquared"))

# Test 1: Algebra is a multiply and subtract

modelOut <- mxRun(model)

objResult <- modelOut[["objective"]]@result
original <- (model[['A']]@values %*% model[['B']]@values - model[['O']]@values) ^ 2
new <- (modelOut[['A']]@values %*% modelOut[['B']]@values - modelOut[['O']]@values) ^ 2
c(objResult, original, new)
