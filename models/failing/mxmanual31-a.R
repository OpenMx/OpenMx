library(OpenMx)

A <- mxMatrix(values = 0.5, nrow = 2, ncol = 1, 
	free = TRUE, name = "A")

D <- mxMatrix(type = "Diag", values = c(0, 0.5), 
	free = c(FALSE, TRUE), nrow = 2, name = "D")

expectedCov <- mxAlgebra(A %*% t(A) + D, "expectedCov")

observedCov <- mxData(matrix(c(1.2, 0.8, 0.8, 1.3)), 'cov', numObs = 150, name = 'observedCov')

objective <- mxMLObjective(covariance = "expectedCov",
	name = "objective")

model <- mxModel(A, D, expectedCov, objective, observedCov)

model <- mxRun(model)

print(model)
