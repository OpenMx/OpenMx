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
	oList <- omxConvertObjectives(flatModel)
	aList <- append(aList, oList)
	bounds <- c()
	state <- c()
	objective <- omxObjectiveIndex(flatModel)
	data <- flatModel@datasets
	output <- .Call("callNPSOL", objective, startVals, 
		bounds, mList, pList, aList, data, state)
	model@output <- output
	model <- omxUpdateModelValues(model, flatModel, output$estimate)
	model <- omxUpdateModelAlgebras(model, flatModel, output$algebras)
	return(model)
}