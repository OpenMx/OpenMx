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

objResult <- mxEvaluate(objective, modelOut)
original <- mxEvaluate((A %*% B - O) ^ 2, model)
new <- mxEvaluate((A %*% B - O) ^ 2, modelOut)
c(objResult, original, new)

omxCheckCloseEnough(mxEvaluate(B[1,1], modelOut), -5, epsilon = 10 ^ -4)
