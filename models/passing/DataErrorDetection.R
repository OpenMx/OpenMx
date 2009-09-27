require(OpenMx)
data <- mxData(type = 'raw', matrix(".", 3, 3, dimnames = list(NULL,c('a','b','c'))))
covariance <- mxMatrix('Symm', 3, 3, values = c(1:6), name = 'cov')
means <- mxMatrix('Full', 1, 3, values = c(1:3), name = 'means')
objective <- mxFIMLObjective('cov', 'means')
model <- mxModel('model', objective, covariance, means, data)
omxCheckError(mxRun(model), paste("The data object", omxQuotes("model.data"),
	"contains an observed matrix that is not of type 'double'"))