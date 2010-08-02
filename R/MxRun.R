#
#   Copyright 2007-2010 The OpenMx Project
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

mxRun <- function(model, ..., intervals = FALSE, silent = FALSE, 
		suppressWarnings = FALSE, unsafe = FALSE,
		checkpoint = FALSE, useSocket = FALSE) {
	if(!silent) cat("Running", model@name, "\n")
	frontendStart <- Sys.time()
	garbageArguments <- list(...)
	if (length(garbageArguments) > 0) {
		stop("mxRun does not accept values for the '...' argument")
	}
	runHelper(model, frontendStart, intervals,
		silent, suppressWarnings, unsafe,
		checkpoint, useSocket)
}

runHelper <- function(model, frontendStart, 
		intervals, silent, suppressWarnings, 
		unsafe, checkpoint, useSocket) {
	omxCheckMatrices(model)
	omxVerifyModel(model)
	dataList <- generateDataList(model)
	dshare <- shareData(model)
	independents <- getAllIndependents(dshare)
	indepTimeStart <- Sys.time()
	independents <- omxLapply(independents, mxRun, 
		intervals = intervals, silent = silent, 
		suppressWarnings = suppressWarnings, unsafe = unsafe,
		checkpoint = checkpoint, useSocket = useSocket)
	indepTimeStop <- Sys.time()
	indepElapsed <- indepTimeStop - indepTimeStart
	if (modelIsHollow(model)) {
		return(processHollowModel(model, independents, 
			dataList, frontendStart, indepElapsed))
	}
	frozen <- lapply(independents, omxFreezeModel)
	model <- omxReplaceModels(model, frozen)
	namespace <- omxGenerateNamespace(model)
	flatModel <- omxFlattenModel(model, namespace)	
	omxCheckNamespace(model, namespace)
	freeFixedValues <- omxCheckVariables(flatModel, namespace)
	flatModel <- convertAlgebras(flatModel, list(startvals=freeFixedValues, 
		values=namespace$values, parameters=namespace$parameters))
	defVars <- generateDefinitionList(flatModel)
	model <- convertDatasets(model, defVars, model@options)
	translation <- translateObjectives(model, namespace, flatModel)
	if (slot(translation, ".estimation")) {
		estimates <- runHelper(translation, frontendStart, 
			intervals, silent, suppressWarnings, 
			unsafe, checkpoint, useSocket)
		# TODO: begin populate the model with free parameter estimates
		# TODO: end populate the model with free parameter estimates
		slot(model, ".postEstimation") <- TRUE 
		model <- translateObjectives(model, namespace, flatModel)
		slot(model, ".postEstimation") <- FALSE 
	} else {
		model <- translation
	}
	# Regenerate the namespace and flatModel
	model <- convertSquareBracketLabels(model)
	namespace <- omxGenerateNamespace(model)
	flatModel <- omxFlattenModel(model, namespace)
	freeFixedValues <- omxCheckVariables(flatModel, namespace)
	oldFlatModel <- flatModel
	flatModel <- constraintsToAlgebras(flatModel)
	flatModel <- convertAlgebras(flatModel, list(startvals=freeFixedValues, 
		values=namespace$values, parameters=namespace$parameters))
	cycleDetection(flatModel)
	flatModel <- populateDefInitialValues(flatModel)
	flatModel <- checkEvaluation(model, flatModel, oldFlatModel)
	parameters <- generateParameterList(flatModel)
	matrices <- generateMatrixList(flatModel)
	algebras <- generateAlgebraList(flatModel)
	startVals <- generateValueList(matrices, parameters)
	defVars <- generateDefinitionList(flatModel)		
	objectives <- convertObjectives(flatModel, model, defVars)
	data <- flatModel@datasets
	algebras <- append(algebras, objectives)
	constraints <- convertConstraints(flatModel)
	intervalList <- generateIntervalList(flatModel, intervals, model@name, parameters)
	communication <- generateCommunicationList(model@name, checkpoint, useSocket, model@options)
	state <- c()
	objective <- getObjectiveIndex(flatModel)
	options <- generateOptionsList(model@options)
	frontendStop <- Sys.time()
	frontendElapsed <- (frontendStop - frontendStart) - indepElapsed
	output <- .Call("callNPSOL", objective, startVals, 
		constraints, matrices, parameters, 
		algebras, data, intervalList, communication, options, state, PACKAGE = "OpenMx")
	backendStop <- Sys.time()
	backendElapsed <- backendStop - frontendStop
	independents <- lapply(independents, undoDataShare, dataList)
	model <- omxReplaceModels(model, independents)
	model <- updateModelMatrices(model, flatModel, output$matrices)
	model <- updateModelAlgebras(model, flatModel, output$algebras)
	model <- undoSquareBracketLabels(model)
	model <- resetDataSortingFlags(model)
	model@output <- processOptimizerOutput(suppressWarnings, flatModel,
		names(matrices), names(algebras),
		names(parameters), names(intervalList), unsafe, output)
	model <- populateRunStateInformation(model, parameters, matrices, 
		objectives, data, flatModel@constraints, independents, defVars)
	frontendStop <- Sys.time()
	frontendElapsed <- frontendElapsed + (frontendStop - backendStop)
	model@output <- calculateTiming(model@output, frontendElapsed,
		backendElapsed, indepElapsed, frontendStop, independents)
	return(model)		
}

