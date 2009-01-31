mxJobRun <- function(model) {
	pList <- omxGenerateParameterList(model)
	mList <- omxGenerateSimpleMatrixList(model)
	aList <- omxGenerateAlgebraList(model)
	startVals <- omxGenerateValueList(model)
	bounds <- c()
	state <- c()
	objective <- omxObjFunConvert(model@objective, model)
	output <- .Call("callNPSOL", objective, startVals, 
		bounds, mList, pList, aList, model@data, state)
	model@output <- output
	model <- omxUpdateModelValues(model, output$estimate)
	return(model)
}