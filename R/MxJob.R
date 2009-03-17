mxJobRun <- function(model) {
	cat("Running", model@name, "\n")
	omxCheckNamespace(model)
	omxCheckMatrices(model)
	dshare <- omxShareData(model)
	independents <- omxGetIndependents(dshare)
	independents <- sfLapply(independents, mxJobRun)
	independents <- lapply(independents, omxFreezeModel)
	depModel <- omxReplaceModels(model, independents)
	flatModel <- omxFlattenModel(depModel)
	data <- omxRemoveDataAliases(flatModel@datasets)
	defNames <- omxGenerateDefinitionNames(data)
	parameters <- omxGenerateParameterList(flatModel, defNames)
	definitions <- omxGenerateDefinitionList(flatModel, defNames)
	matrices <- omxGenerateSimpleMatrixList(flatModel)
	algebras <- omxGenerateAlgebraList(flatModel)
	startVals <- omxGenerateValueList(flatModel, defNames)
	objectives <- omxConvertObjectives(flatModel, definitions)
	algebras <- append(algebras, objectives)
	constraints <- omxConvertConstraints(flatModel)
	state <- c()
	objective <- omxObjectiveIndex(flatModel)
	output <- .Call("callNPSOL", objective, startVals, 
		constraints, matrices, parameters, 
		algebras, data, state)
	model@output <- output
	model <- omxUpdateModelValues(model, 
		flatModel, parameters, output$estimate)
	model <- omxUpdateModelAlgebras(model, 
		flatModel, output$algebras)
	return(model)
}
