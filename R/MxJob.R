setClass(Class = "MxJob",
	representation = representation(
		model = "MxModel",
		objective = "MxObjective",
		data = "MxData"))
	
setMethod("initialize", "MxJob",
	function(.Object, model, objective, data) {
		.Object@model <- model
		.Object@objective <- objective
		.Object@data <- data
		return(.Object)	
	}
)

omxJob <- function(model, objective, data = matrix()) {
	return(new("MxJob", model, objective, data))
}

mxJobRun <- function(job) {
	pList <- omxGenerateParameterList(job@model)
	mList <- omxGenerateSimpleMatrixList(job@model)
	aList <- omxGenerateAlgebraList(job@model)
	startVals <- omxGenerateValueList(job@model)
	bounds <- c()
	state <- c()
	result <- .Call("callNPSOL", job@objective, startVals, 
		bounds, mList, pList, aList, job@data, state)
	return(result)
}



mxJob <- function(model, objective, data = matrix()) {
	return(omxJob(model, objective, data))
}