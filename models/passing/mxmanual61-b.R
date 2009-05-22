require(OpenMx)
W <- mxMatrix("Symm", nrow = 3, ncol = 3, name = "W",
	values = c(1.0, 0.4, 0.3, 0.9, 0.5, 1.1))
H <- mxMatrix("Symm", nrow = 3, ncol = 3, name = "H",
	values = c(1.2, 0.42, 0.3, 1.0, 0.47, 0.9))
M <- mxMatrix("Full", nrow = 3, ncol = 3, name = "M",
	values = c(0.4, 0.05, 0.22, 0.1, 0.3, 0.11, 0.2, 0.12, 0.5))
D <- mxAlgebra(solve(W) %*% M %*% solve(H), name = "D")
model <- mxModel(W, H, M, D)
model <- mxRun(model)
print(model)
