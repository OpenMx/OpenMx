
require(OpenMx)
data(demoOneFactor)
demoOneFactorMiss <- as.matrix(demoOneFactor)[1:5,]
nvar <- ncol(demoOneFactor)
varnames <- colnames(demoOneFactor)



missmat <- matrix(c(2, 4, 4, 4, 4, 4, 5, 1, 1:5, 5), nrow=7, ncol=2)
demoOneFactorMiss[missmat] <- NA

ssModel <- mxModel(name="State Space Missing Debug",
	mxMatrix("Full", 1, 1, FALSE, .3, name="A"),
	mxMatrix("Zero", 1, 1, name="B"),
	mxMatrix("Full", nvar, 1, FALSE, c(.4, .5, .6, .7, .8), name="C", dimnames=list(varnames, "F1")),
	mxMatrix("Zero", nvar, 1, name="D"),
	mxMatrix("Diag", 1, 1, FALSE, 1, name="Q"),
	mxMatrix("Diag", nvar, nvar, FALSE, .2, name="R"),
	mxMatrix("Zero", 1, 1, name="x0"),
	mxMatrix("Diag", 1, 1, FALSE, 1, name="P0"),
	mxData(observed=demoOneFactorMiss, type="raw"),
	imxExpectationStateSpace("A", "B", "C", "D", "Q", "R", "x0", "P0"),
	mxFitFunctionML()
)

ssRun <- mxRun(ssModel)


#------------------------------------------------------------------------------
A <- ssModel$A@values
B <- ssModel$B@values
C <- ssModel$C@values
D <- ssModel$D@values
Q <- ssModel$Q@values
R <- ssModel$R@values
x <- ssModel$x0@values
P <- ssModel$P0@values
u <- matrix(0, 1, 1)

for(i in 1:nrow(demoOneFactorMiss)){
	x <- A %*% x + B %*% u
	print(paste("***************Row***************", i))
	print("Predicted x")
	print(x)
	P <- A %*% P %*% t(A) + Q
	print("Predicted P")
	print(P)
	
	toRemove <- !is.na(demoOneFactorMiss[i,])
	if(sum(!toRemove) < nvar){
		S <- C[toRemove,] %*% P %*% t(C[toRemove,]) + R[toRemove, toRemove]
		print("S")
		print(S)
		S <- solve(S)
		print("S^-1")
		print(S)
		K <- t(P %*% t(C[toRemove,]) %*% S)
		print("K^T")
		print(K)
		print("r")
		y <- demoOneFactorMiss[i,toRemove]
		r <- y - C[toRemove,] %*% x - D[toRemove,] %*% u
		print(r)
		x <- x + t(K)%*%r
		print("Updated x")
		print(x)
		P <- P - t(K) %*% C[toRemove,] %*% P
		print("Updated P")
		print(P)
	}
	else{
		print("S is zero")
		print("S^-1 is zero")
		print("K is zero")
		print("r is zero")
		print("Updated x")
		print(x)
		print("Updated P")
		print(P)
	}
}








