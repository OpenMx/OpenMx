library(OpenMx)

A <- mxMatrix(values = 0.5, nrow = 2, ncol = 1, 
	free = TRUE, name = "A")

D <- mxMatrix(type = "Diag", values = c(0, 0.5), 
	specification = c(0, NA), nrow = 2, free = TRUE, 
	name = "D")

algebra <- mxAlgebra(A %*% t(A) + D, "algebra")

objective <- mxAlgebraObjective("algebra", "objective")

model <- mxModel(A, D, algebra, objective)