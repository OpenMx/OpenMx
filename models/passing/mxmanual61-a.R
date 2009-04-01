library(OpenMx)
A <- mxMatrix("Symm", 
	values = c(1.0, 0.23, 0.34, 2.0, 0.45, 3.0),
	nrow = 3, ncol = 3, name = "A")
B <- mxAlgebra(solve(A), name = "B")
model <- mxModel(A, B)
model <- mxRun(model)
print(model)
