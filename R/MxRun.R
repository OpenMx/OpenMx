#
#   Copyright 2007-2009 The OpenMx Project
#
#   Licensed under the Apache License, Version 2.0 (the "License");
#   you may not use this file except in compliance with the License.
#   You may obtain a copy of the License at
# 
#        http://www.apache.org/licenses/LICENSE-2.0
# 
#   Unless required by applicable law or agreed to in writing, software
#   distributed under the License is distributed on an "AS IS" BASIS,
#   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#   See the License for the specific language governing permissions and
#   limitations under the License.


mxRun <- function(model) {
	cat("Running", model@name, "\n")
	namespace <- omxGenerateNamespace(model)
	omxCheckNamespace(model, namespace)
	omxCheckMatrices(model)
	dshare <- shareData(model)
	independents <- omxGetIndependents(dshare)
	independents <- sfLapply(independents, mxRun)
	independents <- lapply(independents, omxFreezeModel)
	depModel <- omxReplaceModels(model, independents)
	flatModel <- omxFlattenModel(depModel, namespace)
	data <- flatModel@datasets
	omxCheckFreeVariables(flatModel, namespace)
	parameters <- generateParameterList(flatModel)
	definitions <- generateDefinitionList(flatModel)
	matrices <- generateSimpleMatrixList(flatModel)
	algebras <- generateAlgebraList(flatModel)
	startVals <- generateValueList(flatModel, matrices, parameters)
	objectives <- convertObjectives(flatModel, definitions)
	algebras <- generateAlgebraReferences(flatModel, algebras, objectives)
	constraints <- convertConstraints(flatModel)
	state <- c()
	objective <- getObjectiveIndex(flatModel)
	output <- .Call("callNPSOL", objective, startVals, 
		constraints, matrices, parameters, 
		algebras, data, state, PACKAGE = "OpenMx")
	model <- updateModelValues(model, flatModel, parameters, output$estimate)
	model <- updateModelAlgebras(model, flatModel, output$algebras)
	model@output <- computeOptimizationStatistics(flatModel, parameters, output)
	return(model)
}

computeOptimizationStatistics <- function(flatModel, parameters, output) {
	names(output$estimate) <- names(parameters)
	objective <- flatModel@objective
	if (!(is.null(objective) || is(objective, "MxAlgebraObjective"))) {
		output[['AIC']] <- output$minimum + 2 * length(parameters)
#		output[['BIC']] <- output$minimum + length(parameters) * log()
	}
	return(output)
}
