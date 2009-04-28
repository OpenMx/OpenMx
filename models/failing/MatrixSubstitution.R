library(OpenMx)
a <- mxMatrix(values = matrix(1, 1, 1), name = 'a')
b <- mxMatrix(values = matrix(2, 1, 1), name = 'b')
c <- mxMatrix(values = matrix(3, 1, 1), name = 'c')
d <- mxMatrix(values = matrix(4, 1, 1), name = 'd')

e <- mxMatrix(labels = c('a','b','c','d'), 
	nrow = 2, ncol = 2, byrow = TRUE, name = 'e')

model <- mxModel(a, b, c, d, e)

model <- mxRun(model)

expected <- matrix(c(1, 2, 3, 4), nrow = 2, ncol = 2, byrow = TRUE)

cat("Correct: ", all(model[['e']]@values == expected), "\n")
