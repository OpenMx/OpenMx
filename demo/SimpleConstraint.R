library(OpenMx)

O <- mxMatrix("Full", c(3), nrow=1, ncol=1, name="O")
A <- mxMatrix("Full", c(1), nrow=1, ncol=1, name="A")
AlgA <- mxAlgebra(O, name="AlgA")

model <- mxModel()
model <- mxModel(model, O)
model <- mxModel(model, A)
model[["O"]]@specification[1,1] <- "orange"
# Test 1: Algebra is just a matrix.

model <- mxModel(model, AlgA)
model <- mxModel(model, mxAlgebraObjective(algebra = "AlgA"))
model <- mxModel(model, mxConstraint("OBound", "AlgA", "A", ">"))
model <- mxJobRun(model)

model@output
outputA <- model[["AlgA"]]@result
valA <- O[1,1]
diffA <- (valA - outputA[1,1]) / valA
diffA



