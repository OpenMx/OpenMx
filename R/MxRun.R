#
#   Copyright 2007-2013 The OpenMx Project
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

	model <- imxPreprocessModel(model)
	model <- eliminateObjectiveFunctions(model)
	imxCheckMatrices(model)
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
	flatModel <- constraintsToAlgebras(flatModel)
	flatModel <- convertAlgebras(flatModel, convertArguments)
	defVars <- generateDefinitionList(flatModel, list())
	model <- expectationFunctionAddEntities(model, flatModel, labelsData)
	model <- convertDatasets(model, defVars, model@options)
	flatModel@datasets <- collectDatasets(model)
	labelsData <- imxGenerateLabels(model)

	model <- fitFunctionAddEntities(model, flatModel, labelsData)

	if (model@.newobjects) {
		namespace <- imxGenerateNamespace(model)
		flatModel <- imxFlattenModel(model, namespace)
		labelsData <- imxGenerateLabels(model)
	}

	flatModel <- expectationFunctionConvertEntities(flatModel, namespace, labelsData)

	if (model@.newobjects) {
		convertArguments <- imxCheckVariables(flatModel, namespace)
		flatModel <- constraintsToAlgebras(flatModel)
		flatModel <- convertAlgebras(flatModel, convertArguments)
	}

	dependencies <- cycleDetection(flatModel)
	dependencies <- transitiveClosure(flatModel, dependencies)
	flatModel <- populateDefInitialValues(flatModel)
	flatModel <- checkEvaluation(model, flatModel)
	flatModel <- generateFreeVarGroups(flatModel)
	flatModel <- generateParameterList(flatModel, dependencies)
	matrices <- generateMatrixList(flatModel)
	algebras <- generateAlgebraList(flatModel)
	defVars <- generateDefinitionList(flatModel, dependencies)		
	expectations <- convertExpectationFunctions(flatModel, model, labelsData, defVars, dependencies)
	fitfunctions <- convertFitFunctions(flatModel, model, labelsData, defVars, dependencies)
	data <- flatModel@datasets
	numAlgebras <- length(algebras)
	algebras <- append(algebras, fitfunctions)
	constraints <- convertConstraints(flatModel)
	parameters <- flatModel@parameters
	numParam <- length(parameters)
	intervalList <- generateIntervalList(flatModel, intervals, model@name, parameters, labelsData)
	communication <- generateCommunicationList(model@name, checkpoint, useSocket, model@options)

	useOptimizer <- useOptimizer && imxPPML.Check.UseOptimizer(model@options$UsePPML)
	options <- generateOptionsList(model, numParam, constraints, useOptimizer)
	
	compute <- NULL
	computes <- list()
	if (!is.null(model@fitfunction) && is.null(model@compute)) {
		# horrible hack, sorry
		fitNum <- match(flatModel@fitfunction@name, names(flatModel@fitfunctions)) - 1L + numAlgebras
		if (!useOptimizer || numParam == 0) {
			computes <- list(mxComputeOnce(what=fitNum))
		} else {
			if (options[["Calculate Hessian"]] == "No") {
				computes <- list(mxComputeGradientDescent(type="Quasi-Newton",
									  fitfunction=fitNum))
			} else {
				want.se <- options[["Standard Errors"]] == "Yes"
				steps <- list(mxComputeGradientDescent(fitfunction=fitNum, type="Quasi-Newton"),
					      mxComputeEstimatedHessian(fitfunction=fitNum, want.se=want.se))
				computes <- list(mxComputeSequence(steps))
			}
		}
		flatModel@computes <- computes		
		compute <- 0L
	} else {
		if (!is.null(flatModel@compute)) {
			compute <- imxLocateIndex(flatModel, flatModel@compute@name, flatModel@name)
		}
	}

	computes <- convertComputes(flatModel, model)
	
	frontendStop <- Sys.time()
	frontendElapsed <- (frontendStop - frontendStart) - indepElapsed
	if (onlyFrontend) return(model)
	output <- .Call(omxBackend, compute,
			constraints, matrices, flatModel@freeGroupNames, parameters,
			algebras, expectations, computes,
			data, intervalList, communication, options, PACKAGE = "OpenMx")
	backendStop <- Sys.time()
	backendElapsed <- backendStop - frontendStop
	model <- updateModelMatrices(model, flatModel, output$matrices)
	model <- updateModelAlgebras(model, flatModel, output$algebras)
	model <- updateModelExpectations(model, flatModel, output$expectations)
	independents <- lapply(independents, undoDataShare, dataList)
	model <- imxReplaceModels(model, independents)
	model <- resetDataSortingFlags(model)
	model@output <- processOptimizerOutput(suppressWarnings, flatModel,
		names(matrices), names(algebras),
		names(parameters), names(intervalList), output)
	model <- populateRunStateInformation(model, parameters, matrices, 
		fitfunctions, expectations, data, flatModel@constraints, independents, defVars)
	frontendStop <- Sys.time()
	frontendElapsed <- frontendElapsed + (frontendStop - backendStop)
	model@output <- calculateTiming(model@output, frontendElapsed,
		backendElapsed, indepElapsed, frontendStop, independents)
	processErrorConditions(model, unsafe, suppressWarnings)
	return(model)		
}

