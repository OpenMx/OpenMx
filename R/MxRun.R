#
#   Copyright 2007-2014 The OpenMx Project
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

	if (length(intervals) != 1 ||
		typeof(intervals) != "logical" ||
		is.na(intervals)) {
		stop(paste("'intervals' argument",
			"must be TRUE or FALSE in",
			deparse(width.cutoff = 400L, sys.call())), call. = FALSE)
	}

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
	if (!is.null(model@compute)) model@compute <- assignId(model@compute, 1L, '.')
	flatModel <- imxFlattenModel(model, namespace)	
	omxCheckNamespace(model, namespace)
	convertArguments <- imxCheckVariables(flatModel, namespace)
	freeVarGroups <- buildFreeVarGroupList(flatModel)
	flatModel <- constraintsToAlgebras(flatModel)
	flatModel <- convertAlgebras(flatModel, convertArguments)
	defVars <- generateDefinitionList(flatModel, list())
	model <- expectationFunctionAddEntities(model, flatModel, labelsData)
	model <- preprocessDatasets(model, defVars, model@options) # DEPRECATED
	flatModel@datasets <- collectDatasets(model)  # done in imxFlattenModel, but confusingly do it again
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
	flatModel <- generateParameterList(flatModel, dependencies, freeVarGroups)
	matrices <- generateMatrixList(flatModel)
	algebras <- generateAlgebraList(flatModel)
	defVars <- generateDefinitionList(flatModel, dependencies)		
	expectations <- convertExpectationFunctions(flatModel, model, labelsData, defVars, dependencies)
	fitfunctions <- convertFitFunctions(flatModel, model, labelsData, defVars, dependencies)
	data <- convertDatasets(flatModel@datasets, model, flatModel)
	numAlgebras <- length(algebras)
	algebras <- append(algebras, fitfunctions)
	constraints <- convertConstraints(flatModel)
	parameters <- flatModel@parameters
	numParam <- length(parameters)
	intervalList <- generateIntervalList(flatModel, model@name, parameters, labelsData)
	communication <- generateCommunicationList(model@name, checkpoint, useSocket, model@options)

	useOptimizer <- useOptimizer && PPML.Check.UseOptimizer(model@options$UsePPML)
	options <- generateOptionsList(model, numParam, constraints, useOptimizer)
	
	if (!is.null(model@fitfunction) && is.null(model@compute)) {
		# horrible hack, sorry
		compute <- NULL
		fitNum <- paste(model@name, 'fitfunction', sep=".")
		if (!useOptimizer || numParam == 0) {
			compute <- mxComputeOnce(from=fitNum, 'fit', .is.bestfit=TRUE)
		} else {
			steps = list(mxComputeGradientDescent(fitfunction=fitNum))
			if (intervals) {
				steps <- c(steps, mxComputeConfidenceInterval(fitfunction=fitNum))
			}
			if (options[["Calculate Hessian"]] == "Yes") {
				steps <- c(steps, mxComputeNumericDeriv(fitfunction=fitNum))
			}
			if (options[["Standard Errors"]] == "Yes") {
				steps <- c(steps, mxComputeStandardError())
			}
			compute <- mxComputeSequence(c(steps, mxComputeReportDeriv()))
		}
		compute <- assignId(compute, 1L, '.')
		flatModel@compute <- compute
	}

	computes <- convertComputes(flatModel, model)
	
	frontendStop <- Sys.time()
	frontendElapsed <- (frontendStop - frontendStart) - indepElapsed
	if (onlyFrontend) return(model)
	output <- .Call(backend,
			constraints, matrices, parameters,
			algebras, expectations, computes,
			data, intervalList, communication, options, PACKAGE = "OpenMx")
	backendStop <- Sys.time()
	backendElapsed <- backendStop - frontendStop
	model <- updateModelMatrices(model, flatModel, output$matrices)
	model <- updateModelAlgebras(model, flatModel, output$algebras)
	model <- updateModelExpectations(model, flatModel, output$expectations)
	model <- updateModelData(model, flatModel, output$data)
	model@compute <-updateModelCompute(model, output$computes)
	independents <- lapply(independents, undoDataShare, dataList)
	model <- imxReplaceModels(model, independents)
	model@output <- nameOptimizerOutput(suppressWarnings, flatModel,
		names(matrices), names(algebras),
		names(parameters), names(intervalList), output)

	# runstate preserves the pre-backend state of the model
	runstate <- model@runstate
	runstate$parameters <- parameters
	runstate$matrices <- matrices
	runstate$fitfunctions <- fitfunctions
	runstate$expectations <- expectations
	runstate$datalist <- data
	runstate$constraints <- flatModel@constraints
	runstate$independents <- independents
	runstate$defvars <- names(defVars)
	runstate$compute <- computes
	model@runstate <- runstate

	frontendStop <- Sys.time()
	frontendElapsed <- frontendElapsed + (frontendStop - backendStop)
	model@output <- calculateTiming(model@output, frontendElapsed,
		backendElapsed, indepElapsed, frontendStop, independents)
	processErrorConditions(model, unsafe, suppressWarnings)
	return(model)		
}

