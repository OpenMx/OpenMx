library(OpenMx)

O <- mxMatrix("Full", 1, 1, values = 3, name="O")
A <- mxMatrix("Full", 1, 2, values = c(1, 2), name="A")
B <- mxMatrix("Full", 2, 1, values = c(3, 4), free = c(TRUE, FALSE), name="B")
algebra <- mxAlgebra(A %*% B - O, "algebra")
algSquared <- mxAlgebra(algebra %*% algebra, "algebraSquared")

model <- mxModel()
model <- mxModel(model, A, B, O)
model <- mxModel(model, algebra, algSquared)
model <- mxModel(model, mxAlgebraObjective("algebraSquared"))

# Test 1: Algebra is a multiply and subtract

modelOut <- mxRun(model)

objResult <- mxEvaluate(objective, modelOut)
original <- mxEvaluate((A %*% B - O) ^ 2, model)
new <- mxEvaluate((A %*% B - O) ^ 2, modelOut)
c(objResult, original, new)

omxCheckCloseEnough(mxEvaluate(B[1,1], modelOut), -5, epsilon = 10 ^ -4)
