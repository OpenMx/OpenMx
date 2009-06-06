require(OpenMx)

B <- mxMatrix("Full", values=3, free=TRUE, nrow=1, ncol=1, name="B")
A <- mxMatrix("Full", values=1, nrow=1, ncol=1, name="A")
algebra <- mxAlgebra(B, name="algebra")

model <- mxModel()
model <- mxModel(model, B)
model <- mxModel(model, A)

# Test 1: Algebra is just a matrix.
model <- mxModel(model, algebra)
model <- mxModel(model, mxAlgebraObjective("algebra"))
model <- mxModel(model, mxConstraint("B", ">", "A"))
model <- mxRun(model)

omxCheckCloseEnough(model$A[1,1], model$B[1,1], 0.00001)

