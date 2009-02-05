mxJobRun <- function(model) {
	cat("Running", model@name, "\n")
	omxCheckNamespace(model)
	dshare <- omxShareData(model)
	independents <- omxGetIndependents(dshare)
	independents <- lapply(independents, mxJobRun)
	independents <- lapply(independents, omxFreezeModel)
	depModel <- omxReplaceModels(dshare, independents)
	flatModel <- omxFlattenModel(depModel)
	pList <- omxGenerateParameterList(flatModel)
	mList <- omxGenerateSimpleMatrixList(flatModel)
	aList <- omxGenerateAlgebraList(flatModel)
	startVals <- omxGenerateValueList(flatModel)
	bounds <- c()
	state <- c()
	objective <- omxObjFunConvert(flatModel@objective, model)
	output <- .Call("callNPSOL", objective, startVals, 
		bounds, mList, pList, aList, model@data, state)
	model@output <- output
	model <- omxUpdateModelValues(model, model, output$estimate)
	model <- omxUpdateModelObjective(model, output$minimum)
#	model <- omxUpdateModelAlgebras(model, model, output$algebras)
	return(model)
}