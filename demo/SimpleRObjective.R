A <- mxMatrix(nrow = 2, ncol = 2, values = c(1:4), free = TRUE, name = 'A')

squared <- function(x) { x ^ 2 }

objFunction <- function(model, state) {
	values <- model[['A']]@values 
	return(squared(values[1,1] - 4) + squared(values[1,2] - 3) +
		squared(values[2,1] - 2) + squared(values[2,2] - 1))
}
objective <- mxRObjective(objFunction)

model <- mxModel('model', A, objective)

modelOut <- mxRun(model)

omxCheckCloseEnough(modelOut[['A']]@values, c(4, 2, 3, 1), epsilon = 0.001)