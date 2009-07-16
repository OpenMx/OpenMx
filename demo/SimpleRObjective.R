A <- mxMatrix(nrow = 2, ncol = 2, values = c(1:4), free = TRUE, name = 'A')

squared <- function(x) { x ^ 2 }

objFunction <- function(model, state) {
	return(mxEvaluate(squared(A[1,1] - 4) + 
		squared(A[1,2] - 3) +
		squared(A[2,1] - 2) +
		squared(A[2,2] - 1), model))
}
objective <- mxRObjective(objFunction)

model <- mxModel('model', A, objective)

modelOut <- mxRun(model)

omxCheckCloseEnough(mxEvaluate(A, modelOut), 
	rbind(c(4, 3), c(2, 1)), 
	epsilon = 0.001)
