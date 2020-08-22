#
#   Copyright 2007-2019 by the individuals mentioned in the source code history
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
		useOptimizer = TRUE, beginMessage=!silent){

	warnModelCreatedByOldVersion(model)

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
  prohibitDotdotdot(list(...))
	runHelper(model, frontendStart, intervals,
		silent, suppressWarnings, unsafe,
		checkpoint, useSocket, onlyFrontend, useOptimizer, beginMessage)
}

runHelper <- function(model, frontendStart, 
		intervals, silent, suppressWarnings, 
		unsafe, checkpoint, useSocket, onlyFrontend, useOptimizer,
		beginMessage, parentData = NULL) {

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
		onlyFrontend = onlyFrontend, useOptimizer = useOptimizer,
		beginMessage=beginMessage, parentData = model@data)
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
		onlyFrontend = onlyFrontend, beginMessage=beginMessage,
		useOptimizer = useOptimizer)
	indepTimeStop <- Sys.time()
	indepElapsed <- indepTimeStop - indepTimeStart
	if (modelIsHollow(model)) {
		return(processHollowModel(model, independents, 
					  frontendStart, indepElapsed))
	}
	frozen <- lapply(independents, imxFreezeModel)
	model <- imxReplaceModels(model, frozen)
	namespace <- imxGenerateNamespace(model)
	flatModel <- imxFlattenModel(model, namespace, unsafe)
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
		if (is.null(intervals)) intervals <- FALSE
		compute <- omxDefaultComputePlan(modelName=model@name, intervals=intervals,
					useOptimizer=useOptimizer, optionList=options)
		compute@.persist <- FALSE
		model@compute <- compute
	}
	if (is.null(intervals)) intervals <- TRUE
	if (!is.null(model@compute)) model@compute <- assignId(model@compute, 1L, '.')
	flatModelCompute <- safeQualifyNames(model@compute, model@name, namespace)

	omxCheckNamespace(model, namespace)
	convertArguments <- imxCheckVariables(flatModel, namespace)
	flatModel <- constraintsToAlgebras(flatModel)
	flatModel <- eliminateObjectiveFunctions(flatModel)
	flatModel <- convertAlgebras(flatModel, convertArguments)
	defVars <- generateDefinitionList(flatModel, list())
	labelsData <- imxGenerateLabels(model)
	model <- expectationFunctionAddEntities(model, flatModel, labelsData)
	model <- preprocessDatasets(model, defVars, model@options) # DEPRECATED
	flatModel@datasets <- collectDatasets(model, namespace)  # done in imxFlattenModel, but confusingly do it again

	model <- fitFunctionAddEntities(model, flatModel, labelsData)

	if (model@.newobjects) {
		namespace <- imxGenerateNamespace(model)
		flatModel <- imxFlattenModel(model, namespace, unsafe)
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

	intervalList <- NULL
	if (intervals) {
		intervalList <- generateIntervalList(flatModel, model@name, parameters, labelsData)
	}
	communication <- generateCommunicationList(model, checkpoint, useSocket, model@options)

	useOptimizer <- useOptimizer && PPML.Check.UseOptimizer(model@options$UsePPML)
	options <- limitMajorIterations(options, numParam, length(constraints))
	computes <- convertComputes(flatModel, model)
	
	frontendStop <- Sys.time()
	frontendElapsed <- (frontendStop - frontendStart) - indepElapsed
	if(beginMessage) message("Running ", model@name, " with ", numParam, " parameter",
			    ifelse(numParam==1, "", "s"))
	if (onlyFrontend) return(model)

	output <- .Call(backend,
			constraints, matrices, parameters,
			algebras, expectations, computes,
			data, intervalList, communication, options, defVars,
			silent || !interactive(), PACKAGE = "OpenMx")
	backendStop <- Sys.time()
	backendElapsed <- backendStop - frontendStop
	if (is.null(output$error)) {
		model <- updateModelMatrices(model, flatModel, output$matrices)
		model <- updateModelAlgebras(model, flatModel, output$algebras)
		model <- updateModelExpectations(model, flatModel, output$expectations)
		model <- updateModelExpectationDims(model, expectations)
		model <- updateModelData(model, flatModel, output$data)
	}
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
	if (is.na(model@output$status$code) ||
	    (!is.na(model@output$status$code) && model@output$status$code < 5)) {
		if (!is.null(model@output[['infoDefinite']]) &&
		    !is.na(model@output[['infoDefinite']]) && !model@output[['infoDefinite']]) {
			model@output$status$code <- 5   # INFORM_NOT_CONVEX
		}
	}
	mroe <- model@output[['maxRelativeOrdinalError']]
	if (!is.null(mroe) && mroe > options[["mvnRelEps"]]) {
		warning(paste("Polite note: Model finished with a larger ordinal error than we typically expect.\n",
			"This may be fine, but you may wish to re-run the model using\n",
			"`mxTryHardOrdinal()` in place of `mxRun()` to try for a better fit.\n",
			"Expert version: model$output[['maxRelativeOrdinalError']] is \n",
			"larger than the mvnRelEps value of ", options[["mvnRelEps"]], ".\n",
			"If this is expected for your model, you might wish to increase `mvnRelEps`, e.g:\n",
			"mxOption(NULL, 'mvnRelEps', value= mxOption(NULL, 'mvnRelEps')*5)\n",
			"see `?mxOptions`" ))
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

imxReportProgress <- function(info, eraseLen) {
	origLen = nchar(info)
	if (origLen < eraseLen) {
		info <- paste0(info, paste0(rep(' ', eraseLen - nchar(info)), collapse=""))
	}
	cat(paste0("\r", info))
	if (origLen == 0) cat("\r")
}

enumerateDatasets <- function(model) {
	datasets <- c()
	if (!is.null(model@data)) datasets <- c(datasets, model@name)
	if (length(model@submodels)) for (mx in 1:length(model@submodels)) {
    got <- enumerateDatasets(model@submodels[[mx]])
    if (length(got)) datasets <- c(datasets, got)
	}
	return(datasets)
}

as.statusCode <- function(code) {
	lev <- c("OK", "OK/green",
		"infeasible linear constraint",
		"infeasible non-linear constraint",
		"iteration limit/blue",
		"not convex/red",
		"nonzero gradient/red",
		"bad deriv",
		"internal error",
		"infeasible start")
	if (is(code, 'ordered')) {
		if (all(levels(code) == lev)) return(code)
		code <- as.character(code)
	}
	if (is.numeric(code) || is.null(code)) {
		mxFactor(code, levels=c(0:7,9:10), labels=lev)
	} else if (is.character(code) || all(is.na(code))) {
		mxFactor(code, lev)
	} else {
		stop(paste("Don't know how to convert type", typeof(code),
			   "into a status code"))
	}
}

mxBootstrap <- function(model, replications=200, ...,
                        data=NULL, plan=NULL, verbose=0L,
                        parallel=TRUE, only=as.integer(NA),
			OK=mxOption(model, "Status OK"), checkHess=FALSE) {
  warnModelCreatedByOldVersion(model)
    if (!missing(plan)) {
	    stop("The 'plan' argument is deprecated. Use mxModel(model, plan) and then mxBootstrap")
    }
    if (!is(model$compute, "MxComputeBootstrap")) {
	    if (missing(checkHess)) checkHess <- as.logical(NA)
	    model <- ProcessCheckHess(model, checkHess)
    if (missing(data)) {
      data <- enumerateDatasets(model)
    }
	    # wrap with mxComputeTryCatch ?
    plan <- mxComputeBootstrap(data, model@compute)
  } else {
    plan <- model$compute
  }

  plan$replications <- as.integer(replications)
  plan$verbose <- as.integer(verbose)
  plan$parallel <- as.logical(parallel)
  plan$only <- as.integer(only)
  plan$OK <- as.statusCode(OK)
  
  model <- mxModel(model, plan)
  mxRun(model, suppressWarnings=TRUE)
}

omxGetBootstrapReplications <- function(model) {
   if(!is(model, "MxModel")) {
      stop("'model' argument must be a MxModel object")
   }
  if (is.null(model$compute) || !is(model$compute, "MxComputeBootstrap")) {
	  stop(paste("Compute plan", class(model$compute), "found in model",
		     omxQuotes(model$name),
		     "instead of MxComputeBootstrap. You need to run this model through mxBootstrap first."))
  }
   assertModelFreshlyRun(model)
  cb <- model@compute
  if (is.null(cb@output$raw)) {
	  stop(paste("No bootstrap data found. Please run this model",
		     "through mxBootstrap again."))
  }
  if (!is.na(cb@only)) {
	  stop(paste("Detected mxBootstrap's only= option. Please mxBootstrap",
		     "this model without using only="))
  }
  if (cb@output$numParam != length(coef(model))) {
	  stop(paste("Model", omxQuotes(model), "has", length(coef(model)),
		     "parameters but bootstrap data has", cb@output$numParam,
		     "parameters. Please mxBootstrap this model again."))
  }
  raw <- cb@output$raw
  mask <- raw[,'statusCode'] %in% cb@OK
  bootData <- raw[mask, 3:(length(coef(model))+2), drop=FALSE]
  if (sum(mask) < 3) {
	  stop(paste("Fewer than 3 replications are available.",
		     "Use mxBootstrap to increase the number of replications."))
  }
   if (sum(mask) < .95*length(mask)) {
	   pct <- round(100*sum(mask) / length(mask))
	   warning(paste0("Only ",pct,"% of the bootstrap replications ",
			  "converged acceptably. Accuracy is much less than the ", nrow(raw),
			  " replications requested. Examine table(model$compute$output$raw$statusCode)"), call.=FALSE)
   }
   bootData
}

omxBootstrapCov <- function(model) cov(omxGetBootstrapReplications(model))
