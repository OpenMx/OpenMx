#
#   Copyright 2007-2012 The OpenMx Project
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
		checkpoint = FALSE, useSocket = FALSE, onlyFrontend = FALSE, 
		useOptimizer = TRUE){
	if(!silent) cat("Running", model@name, "\n")
	frontendStart <- Sys.time()
	garbageArguments <- list(...)
	if (length(garbageArguments) > 0) {
		stop("mxRun does not accept values for the '...' argument")
	}
	runHelper(model, frontendStart, intervals,
		silent, suppressWarnings, unsafe,
		checkpoint, useSocket, onlyFrontend, useOptimizer)
}

runHelper <- function(model, frontendStart, 
		intervals, silent, suppressWarnings, 
		unsafe, checkpoint, useSocket, onlyFrontend, useOptimizer, parentData = NULL) {
	omxCheckMatrices(model)
	imxVerifyModel(model)
	model <- processParentData(model, parentData)
	if (modelIsHollow(model)) {
		independents <- getAllIndependents(model)
		indepTimeStart <- Sys.time()
	    independents <- omxLapply(independents, runHelper,
		frontendStart = frontendStart, 
		intervals = intervals, silent = silent, 
		suppressWarnings = suppressWarnings, unsafe = unsafe,
		checkpoint = checkpoint, useSocket = useSocket,
		onlyFrontend = onlyFrontend, useOptimizer = useOptimizer, parentData = model@data)
		indepTimeStop <- Sys.time()
		indepElapsed <- indepTimeStop - indepTimeStart
		return(processHollowModel(model, independents, 
			frontendStart, indepElapsed))
	}
	dataList <- generateDataList(model)
	dshare <- shareData(model)
	independents <- getAllIndependents(dshare)
	indepTimeStart <- Sys.time()
	independents <- omxLapply(independents, mxRun, 
		intervals = intervals, silent = silent, 
		suppressWarnings = suppressWarnings, unsafe = unsafe,
		checkpoint = checkpoint, useSocket = useSocket,
		onlyFrontend = onlyFrontend, useOptimizer = useOptimizer)
	indepTimeStop <- Sys.time()
	indepElapsed <- indepTimeStop - indepTimeStart
	if (modelIsHollow(model)) {
		return(processHollowModel(model, independents, 
			dataList, frontendStart, indepElapsed))
	}
	frozen <- lapply(independents, imxFreezeModel)
	model <- imxReplaceModels(model, frozen)
	namespace <- imxGenerateNamespace(model)
	flatModel <- imxFlattenModel(model, namespace)	
	omxCheckNamespace(model, namespace)
	convertArguments <- imxCheckVariables(flatModel, namespace)
	flatModel <- convertAlgebras(flatModel, convertArguments)
	defVars <- generateDefinitionList(flatModel)
	model <- convertDatasets(model, defVars, model@options)
	flatModel@datasets <- collectDatasets(model)	
	labelsData <- imxGenerateLabels(model)
	triple <- translateObjectives(model, namespace, labelsData, flatModel)
	model <- triple[[1]]
	namespace <- triple[[2]]
	flatModel <- triple[[3]]
	labelsData <- imxGenerateLabels(model)
	convertArguments <- imxCheckVariables(flatModel, namespace)
	flatModel <- constraintsToAlgebras(flatModel)
	flatModel <- convertAlgebras(flatModel, convertArguments)
	cycleDetection(flatModel)
	flatModel <- populateDefInitialValues(flatModel)
	flatModel <- checkEvaluation(model, flatModel)
	parameters <- generateParameterList(flatModel)
	matrices <- generateMatrixList(flatModel)
	algebras <- generateAlgebraList(flatModel)
	startVals <- generateValueList(matrices, parameters)
	defVars <- generateDefinitionList(flatModel)		
	objectives <- convertObjectives(flatModel, model, labelsData, defVars)
	data <- flatModel@datasets
	algebras <- append(algebras, objectives)
	constraints <- convertConstraints(flatModel)
	intervalList <- generateIntervalList(flatModel, intervals, model@name, parameters)
	communication <- generateCommunicationList(model@name, checkpoint, useSocket, model@options)
	state <- c()
	objective <- getObjectiveIndex(flatModel)
	options <- generateOptionsList(model, parameters, constraints, useOptimizer)
	frontendStop <- Sys.time()
	frontendElapsed <- (frontendStop - frontendStart) - indepElapsed
	if (onlyFrontend) return(model)
	output <- .Call("callNPSOL", objective, startVals, 
		constraints, matrices, parameters, 
		algebras, data, intervalList, communication, options, state, PACKAGE = "OpenMx")
	backendStop <- Sys.time()
	backendElapsed <- backendStop - frontendStop
	model <- updateModelMatrices(model, flatModel, output$matrices)
	model <- updateModelAlgebras(model, flatModel, output$algebras)
	independents <- lapply(independents, undoDataShare, dataList)	
	model <- imxReplaceModels(model, independents)
	model <- resetDataSortingFlags(model)
	model@output <- processOptimizerOutput(suppressWarnings, flatModel,
		names(matrices), names(algebras),
		names(parameters), names(intervalList), output)
	model <- populateRunStateInformation(model, parameters, matrices, 
		objectives, data, flatModel@constraints, independents, defVars)
	frontendStop <- Sys.time()
	frontendElapsed <- frontendElapsed + (frontendStop - backendStop)
	model@output <- calculateTiming(model@output, frontendElapsed,
		backendElapsed, indepElapsed, frontendStop, independents)
	processErrorConditions(model, unsafe, suppressWarnings)
	return(model)		
}

