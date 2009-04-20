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
	omxCheckNamespace(model)
	omxCheckMatrices(model)
	dshare <- omxShareData(model)
	independents <- omxGetIndependents(dshare)
	independents <- sfLapply(independents, mxRun)
	independents <- lapply(independents, omxFreezeModel)
	depModel <- omxReplaceModels(model, independents)
	flatModel <- omxFlattenModel(depModel)
	data <- flatModel@datasets
	defNames <- omxGenerateDefinitionNames(data)
	omxCheckFreeVariables(flatModel, defNames)
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
		algebras, data, state, PACKAGE="OpenMx")
	model <- omxUpdateModelValues(model, 
		flatModel, parameters, output$estimate)
	model <- omxUpdateModelAlgebras(model, 
		flatModel, output$algebras)
	model@output <- omxComputeOptimizationStatistics(flatModel, parameters, output)
	return(model)
}

omxComputeOptimizationStatistics <- function(flatModel, parameters, output) {
	objective <- flatModel@objective
	if(!(is.null(objective) || is(objective,"MxAlgebraObjective"))) {
		output[['AIC']] <- output$minimum + 2 * length(parameters)
#		output[['BIC']] <- output$minimum + length(parameters) * log()
	}
	return(output)
}
