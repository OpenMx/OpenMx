#
#   Copyright 2007-2017 The OpenMx Project
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

mxRun <- function(model, ..., intervals=NULL, silent = FALSE, 
		suppressWarnings = FALSE, unsafe = FALSE,
		checkpoint = FALSE, useSocket = FALSE, onlyFrontend = FALSE, 
		useOptimizer = TRUE){

	if (.hasSlot(model, '.version')) {
		mV <- package_version(model@.version)
		curV <- packageVersion('OpenMx')
		if (curV$major != mV$major ||
		    curV$minor != mV$minor) {
			warning(paste0("You are using OpenMx version ", curV,
				       " with a model created by OpenMx version ",
				       mV, ". This may work fine (fingers crossed), but if you run into ",
				       "trouble then please recreate your model with the ",
				       "current version of OpenMx."))
		}
	} else {
		curV <- packageVersion('OpenMx')
		warning(paste0("You are using OpenMx version ", curV,
			       " with a model created by an old version of OpenMx. ",
			       "This may work fine (fingers crossed), but if you run into ",
			       "trouble then please recreate your model with the ",
			       "current version of OpenMx."))
	}

	if (is.null(intervals)) {
		# OK
	} else if (length(intervals) != 1 ||
		typeof(intervals) != "logical" ||
		is.na(intervals)) {
		stop(paste("'intervals' argument",
			"must be TRUE or FALSE in",
			deparse(width.cutoff = 400L, sys.call())), call. = FALSE)
	}

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

	Rcpp::Module  # ensure Rcpp is loaded
	model <- imxPreprocessModel(model)
	model <- eliminateObjectiveFunctions(model)
	model <- zapExtraneousMatrices(model)
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
					  frontendStart, indepElapsed))
	}
	frozen <- lapply(independents, imxFreezeModel)
	model <- imxReplaceModels(model, frozen)
	namespace <- imxGenerateNamespace(model)
	flatModel <- imxFlattenModel(model, namespace)	
	options <- generateOptionsList(model, length(flatModel@constraints), useOptimizer)
	options[['intervals']] <- intervals

	if (!is.null(model@compute) && (!.hasSlot(model@compute, '.persist') || !model@compute@.persist)) {
		model@compute <- NULL
	}
	if (!is.null(model@expectation) && is.null(model@fitfunction) && is.null(model@compute)) {
		# The purpose of this check is to prevent analysts new to OpenMx
		# from running nonsensical models.
		stop(paste(model@name, " has expectation ", class(model@expectation),
			   ", but there is no fitfunction given, and no default.\n",
			   "To fix, see, e.g. help(mxFitFunctionML) for an example fit function, and how these pair with the expectation", sep = ""))
	}

	defaultComputePlan <- (is.null(model@compute) || is(model@compute, 'MxComputeDefault'))
	if (!useOptimizer && !defaultComputePlan) {
		warning("mxRun(..., useOptimizer=FALSE) ignored due to custom compute plan")
	}
	if (!is.null(model@fitfunction) && defaultComputePlan) {
		compute <- omxDefaultComputePlan(modelName=model@name, intervals=(length(intervals) && intervals),
					useOptimizer=useOptimizer, optionList=options)
		compute@.persist <- FALSE
		model@compute <- compute
	}
	if (!is.null(model@compute)) model@compute <- assignId(model@compute, 1L, '.')
	flatModelCompute <- safeQualifyNames(model@compute, model@name, namespace)

	omxCheckNamespace(model, namespace)
	convertArguments <- imxCheckVariables(flatModel, namespace)
	flatModel <- constraintsToAlgebras(flatModel)
	flatModel <- eliminateObjectiveFunctions(flatModel)
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
		flatModel <- eliminateObjectiveFunctions(flatModel)
		flatModel <- convertAlgebras(flatModel, convertArguments)
	}

	dependencies <- cycleDetection(flatModel)
	dependencies <- transitiveClosure(flatModel, dependencies)
	flatModel <- populateDefInitialValues(flatModel)
	flatModel <- checkEvaluation(model, flatModel)
	flatModel@compute <- flatModelCompute
	freeVarGroups <- buildFreeVarGroupList(flatModel)
	flatModel <- generateParameterList(flatModel, dependencies, freeVarGroups)
	matrices <- generateMatrixList(flatModel)
	algebras <- generateAlgebraList(flatModel)
	if (length(defVars)) {
		# We're only going to find them if we found them the first time
		defVars <- generateDefinitionList(flatModel, dependencies)
	}
	expectations <- convertExpectationFunctions(flatModel, model, labelsData, dependencies)

	if (length(expectations)) {
		prec <- lapply(expectations, function(x){getPrecisionPerExpectation(x,options)})

		functionPrecision <- Reduce(max, c(imxAutoOptionValue("Function precision",options),
						   sapply(prec, function(x) x[['functionPrecision']])))
		options[['Function precision']] <- as.character(functionPrecision)

		if (defaultComputePlan && is(model@compute, "MxComputeSequence")) {
			iterations <- ifelse(
				is.na(suppressWarnings(as.numeric(options[["Gradient iterations"]]))),
				Reduce(min, c(4L, sapply(prec, function(x) x[['iterations']]))),
				as.integer(options[["Gradient iterations"]]))
			stepSize <- ifelse(
				is.na(suppressWarnings(as.numeric(options[["Gradient step size"]]))),
				Reduce(max, c(sqrt(.Machine$double.eps),sapply(prec, function(x) x[['stepSize']]))),
				as.numeric(options[["Gradient step size"]]))
			model <- adjustDefaultNumericDeriv(model, iterations, stepSize)
			flatModel <- adjustDefaultNumericDeriv(flatModel, iterations, stepSize)
		}
	}
	else{options[["Function precision"]] <- as.character(imxAutoOptionValue("Function precision",options))}

	fitfunctions <- convertFitFunctions(flatModel, model, labelsData, dependencies)
	data <- convertDatasets(flatModel@datasets, model, flatModel)
	numAlgebras <- length(algebras)
	algebras <- append(algebras, fitfunctions)
	constraints <- convertConstraints(flatModel)
	parameters <- flatModel@parameters
	numParam <- length(parameters)
	if (numParam == 0 && defaultComputePlan && !is.null(model@fitfunction)) {
		compute <- mxComputeSequence(list(CO=mxComputeOnce(from=paste(model@name, 'fitfunction', sep="."),
								   'fit', .is.bestfit=TRUE),
						  RE=mxComputeReportExpectation()))
		compute@.persist <- FALSE
		compute <- assignId(compute, 1L, '.')
		model@compute <- compute
		flatModel@compute <- compute
	}

	intervalList <- generateIntervalList(flatModel, model@name, parameters, labelsData)
	communication <- generateCommunicationList(model, checkpoint, useSocket, model@options)

	useOptimizer <- useOptimizer && PPML.Check.UseOptimizer(model@options$UsePPML)
	options <- limitMajorIterations(options, numParam, length(constraints))
	computes <- convertComputes(flatModel, model)
	
	frontendStop <- Sys.time()
	frontendElapsed <- (frontendStop - frontendStart) - indepElapsed
	if(!silent) message("Running ", model@name, " with ", numParam, " parameter",
			    ifelse(numParam==1, "", "s"))
	if (onlyFrontend) return(model)

	output <- .Call(backend,
			constraints, matrices, parameters,
			algebras, expectations, computes,
			data, intervalList, communication, options, defVars,
			silent || !interactive(), PACKAGE = "OpenMx")
	backendStop <- Sys.time()
	backendElapsed <- backendStop - frontendStop
	model <- updateModelMatrices(model, flatModel, output$matrices)
	model <- updateModelAlgebras(model, flatModel, output$algebras)
	model <- updateModelExpectations(model, flatModel, output$expectations)
	model <- updateModelExpectationDims(model, expectations)
	model <- updateModelData(model, flatModel, output$data)
	model@compute <-updateModelCompute(model, output$computes)
	output[['computes']] <- NULL
	if (!is.null(output[['bounds']])) {
		model <- omxSetParameters(model, names(parameters),
					  lbound=output[['bounds']][['l']],
					  ubound=output[['bounds']][['u']])
		output[['bounds']] <- NULL
	}
	independents <- lapply(independents, undoDataShare, dataList)
	model <- imxReplaceModels(model, independents)
	model@output <- nameOptimizerOutput(suppressWarnings, flatModel,
		names(matrices), names(algebras),
		names(parameters), names(constraints), model@compute, output)
	
	theFitUnits <- model$output$fitUnits
	if(length(theFitUnits) > 0 && theFitUnits %in% "r'Wr" ){
		if(options[["Standard Errors"]] == "Yes"){
			wlsSEs <- imxWlsStandardErrors(model)
			model@output$standardErrors <- wlsSEs$SE
			model@output$hessian <- 2*solve(wlsSEs$Cov) #puts in same units as m2ll Hessian
			wlsChi <- imxWlsChiSquare(model, J=wlsSEs$Jac)
		} else {
			wlsChi <- list(Chi=NA, ChiDoF=NA)
		}
		model@output$chi <- wlsChi$Chi
		model@output$chiDoF <- wlsChi$ChiDoF
	}
	if (is.na(model@output$status$code) ||
	    (!is.na(model@output$status$code) && model@output$status$code < 5)) {
		if (!is.null(model@output[['infoDefinite']]) &&
		    !is.na(model@output[['infoDefinite']]) && !model@output[['infoDefinite']]) {
			model@output$status$code <- 5   # INFORM_NOT_CONVEX
		}
	}

	# Currently runstate preserves the pre-backend state of the model.
	# Eventually this needs to capture the post-backend state,
	# but we need tests in place for summary output to ensure that we
	# don't cause regressions.

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

	model <- clearModifiedSinceRunRecursive(model)

	return(model)		
}

updateModelExpectationDims <- function(model, expectations){
	expectationNames <- names(expectations)
	for(aname in expectationNames){
		if(!is.null(model[[aname]])){
			if (.hasSlot(expectations[[aname]], 'dims')) {
				model[[aname]]@.runDims <- expectations[[aname]]@dims
			}
		}
	}
	return(model)
}

#' Report backend progress
#'
#' Prints a show status string to the console without emitting a
#' newline.
#'
#' @param info the character string to print
#' @param eraseLen the number of characters to erase
imxReportProgress <- function(info, eraseLen) {
	origLen = nchar(info)
	if (origLen < eraseLen) {
		info <- paste0(info, paste0(rep(' ', eraseLen - nchar(info)), collapse=""))
	}
	cat(paste0("\r", info))
	if (origLen == 0) cat("\r")
}

enumerateExpectations <- function(model) {
	expectations <- c()
	if (!is.null(model@expectation)) expectations <- c(expectations, model@name)
	if (length(model@submodels)) {
		expectations <- c(expectations, sapply(model@submodels, enumerateExpectations))
	}
	return(expectations)
}

mxBootstrap <- function(model, replications=200, quantile=c(.25,.75), ...,
                        expectation=NULL, plan=NULL, verbose=0L,
                        parallel=TRUE, seed=42L, only=as.integer(NA)) {
  if (!is(model$compute, "MxComputeBootstrap")) {
    if (!is.null(plan)) {
      plan <- model$compute
    }
    if (missing(expectation)) {
      expectation <- enumerateExpectations(model)
    }
    plan <- mxComputeBootstrap(expectation, plan)
  } else {
    if (!missing(plan)) stop(paste("Model", omxQuotes(model@name), "already has",
                                   "a", omxQuotes(class(model$class)), "plan"))
    plan <- model$compute
  }

  plan$replications <- as.integer(replications)
  plan$quantile <- quantile
  plan$seed <- as.integer(seed)
  plan$verbose <- as.integer(verbose)
  plan$parallel <- as.logical(parallel)
  plan$seed <- as.integer(seed)
  plan$only <- as.integer(only)
  
  model <- mxModel(model, plan)
  mxRun(model)
}
