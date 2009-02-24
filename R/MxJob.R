mxJobRun <- function(model) {
	cat("Running", model@name, "\n")
	omxCheckNamespace(model)
	dshare <- omxShareData(model)
	independents <- omxGetIndependents(dshare)
	independents <- lapply(independents, mxJobRun)
	independents <- lapply(independents, omxFreezeModel)
	depModel <- omxReplaceModels(model, independents)
	flatModel <- omxFlattenModel(depModel)
	parameters <- omxGenerateParameterList(flatModel)
	matrices <- omxGenerateSimpleMatrixList(flatModel)
	algebras <- omxGenerateAlgebraList(flatModel)
	startVals <- omxGenerateValueList(flatModel)
	objectives <- omxConvertObjectives(flatModel)
	algebras <- append(algebras, objectives)
	constraints <- omxConvertConstraints(flatModel)
	state <- c()
	objective <- omxObjectiveIndex(flatModel)
	data <- omxRemoveDataAliases(flatModel@datasets)
	output <- .Call("callNPSOL", objective, startVals, 
		constraints, matrices, parameters, 
		algebras, data, state)
	model@output <- output
	model <- omxUpdateModelValues(model, 
		flatModel, output$estimate)
	model <- omxUpdateModelAlgebras(model, 
		flatModel, output$algebras)
	return(model)
}