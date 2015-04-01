#
#   Copyright 2007-2015 The OpenMx Project
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

	if(!silent) message("Running ", model@name)
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
	if (!is.null(model@compute)) model@compute <- assignId(model@compute, 1L, '.')
	flatModel <- imxFlattenModel(model, namespace)	
	omxCheckNamespace(model, namespace)
	convertArguments <- imxCheckVariables(flatModel, namespace)
	freeVarGroups <- buildFreeVarGroupList(flatModel)
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
	flatModel <- generateParameterList(flatModel, dependencies, freeVarGroups)
	matrices <- generateMatrixList(flatModel)
	algebras <- generateAlgebraList(flatModel)
	if (length(defVars)) {
		# We're only going to find them if we found them the first time
		defVars <- generateDefinitionList(flatModel, dependencies)
	}
	expectations <- convertExpectationFunctions(flatModel, model, labelsData, defVars, dependencies)
	fitfunctions <- convertFitFunctions(flatModel, model, labelsData, defVars, dependencies)
	data <- convertDatasets(flatModel@datasets, model, flatModel)
	numAlgebras <- length(algebras)
	algebras <- append(algebras, fitfunctions)
	constraints <- convertConstraints(flatModel)
	parameters <- flatModel@parameters
	numParam <- length(parameters)
	intervalList <- generateIntervalList(flatModel, model@name, parameters, labelsData)
	communication <- generateCommunicationList(model, checkpoint, useSocket, model@options)

	useOptimizer <- useOptimizer && PPML.Check.UseOptimizer(model@options$UsePPML)
	options <- generateOptionsList(model, numParam, constraints, useOptimizer)
	
	defaultCompute <- NULL
	if (!is.null(model@expectation) && is.null(model@fitfunction) && is.null(model@compute)) {
		# The purpose of this check is to prevent analysts new to OpenMx
		# from running nonsensical models.
		stop(paste(model@name, " has expectation ", class(model@expectation),
			   ", but there is no fitfunction given, and no default.\n",
			   "To fix, see, e.g. help(mxFitFunctionML) for an example fit function, and how these pair with the expectation", sep = ""))
	}
	if (!is.null(model@fitfunction) && is.null(model@compute)) {
		# horrible hack, sorry
		compute <- NULL
		fitNum <- paste(model@name, 'fitfunction', sep=".")
		if (!useOptimizer || numParam == 0) {
			compute <- mxComputeOnce(from=fitNum, 'fit', .is.bestfit=TRUE)
		} else {
			steps = list(mxComputeGradientDescent(fitfunction=fitNum))
			if (intervals) {
				steps <- c(steps, mxComputeConfidenceInterval(
				    mxComputeGradientDescent(fitfunction=fitNum, nudgeZeroStarts=FALSE)))
			}
			if (options[["Calculate Hessian"]] == "Yes") {
				steps <- c(steps, mxComputeNumericDeriv(fitfunction=fitNum))
			}
			if (options[["Standard Errors"]] == "Yes") {
				steps <- c(steps, mxComputeStandardError(), mxComputeHessianQuality())
			}
			compute <- mxComputeSequence(c(steps, mxComputeReportDeriv()))
		}
		compute <- assignId(compute, 1L, '.')
		flatModel@compute <- compute
		defaultCompute <- compute
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
	model <- updateModelExpectationDims(model, expectations)
	model <- updateModelData(model, flatModel, output$data)
	model@compute <-updateModelCompute(model, output$computes)
	independents <- lapply(independents, undoDataShare, dataList)
	model <- imxReplaceModels(model, independents)
	model@output <- nameOptimizerOutput(suppressWarnings, flatModel,
		names(matrices), names(algebras),
		names(parameters), output)
	if(length(model$data) > 0 && model$data$type=="acov"){
		wlsSEs <- imxWlsStandardErrors(model)
		model@output$standardErrors <- wlsSEs$SE
		model@output$hessian <- wlsSEs$Cov
		model@output$calculatedHessian <- wlsSEs$Cov
		wlsChi <- imxWlsChiSquare(model, J=wlsSEs$Jac)
		model@output$chi <- wlsChi$Chi
		model@output$chiDoF <- wlsChi$ChiDoF
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
	if (!is.null(defaultCompute)) {
		runstate$compute <- defaultCompute  # more hack
	} else {
		runstate$compute <- computes
	}
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
			model[[aname]]@.runDims <- expectations[[aname]]@dims
		}
	}
	return(model)
}


# mxTryHard is Wrapper to mxRun that makes multiple attempts to reach an acceptable solution.
# possible TODO's:
#   *Randomly disturb start values very close to zero by adding a random number to them, instead of just multiplying 
#     by a random number.
#   *Edit function so that it does not go through all attempts if model supplied is bad (unidentifed, etc.)?
#   *Stop further tries if fit function value is getting worse, or is improving by less than some amount (Mike Neale's
#     deltas).

mxTryHard<-function (model, extraTries = 10, greenOK = FALSE, loc = 1, 
  scale = 0.25, checkHess = TRUE, fit2beat = Inf, paste = TRUE,
  iterationSummary=FALSE, bestInitsOutput=TRUE, showInits=FALSE,
  ...) 
{
  
  stopflag <- FALSE
  numdone <- 0
  bestfitsofar<-Inf
  inits<-omxGetParameters(model)
  while (!stopflag) {
    if(iterationSummary==TRUE) message(paste0('\nBegin fit attempt ', numdone+1, ' of at maximum ', extraTries +1, ' tries'))
    params <- omxGetParameters(model)
    if(showInits==TRUE) {
      message('Starting values:  ')
      message(paste0(names(params),' : ', params,'\n'))
    }
    numdone <- numdone + 1
    
    fit <- suppressWarnings(try(mxRun(model, suppressWarnings = T, unsafe=T, silent=T,intervals=FALSE)))
    if (class(fit) == "try-error" || fit$output$status$status== -1) {
      newparams<-omxGetParameters(model) #get recent fit
      if(exists('bestfit')) newparams<-bestfit.params #if bestfit exists use this instead
      if(numdone %% 4 == 0) newparams<-inits #sometimes, use initial start values instead
      #       if(numdone %% 5 == 0) { #sometimes, switch optimizers
      #         if(mxOption(NULL, "Default optimizer")=='CSOLNP') newoptimizer<-'NPSOL'
      #         if(mxOption(NULL, "Default optimizer")=='NPSOL') newoptimizer<-'CSOLNP'
      #         message(paste0('Switching to ',newoptimizer,' optimizer for OpenMx temporarily')) 
      #         mxOption(NULL, "Default optimizer", newoptimizer)
      #       }
      model <- omxSetParameters(model, labels = names(newparams), 
        values = newparams * runif(length(params),loc-scale,loc+scale))  #set to multiply bestfit.params instead of params
    }
    else { #if fit was not an error
      if (fit$output$minimum <= bestfitsofar) {
        bestfit <- fit
        bestfit.params <- omxGetParameters(bestfit)
      }
      
      if(fit$output$minimum < bestfitsofar) bestfitsofar <- fit$output$minimum
      
      if (length(fit$output$calculatedHessian) == 0) {
        checkHess <- FALSE
      }
      if (checkHess) {
        if (sum(is.na(fit$output$calculatedHessian)) > 
            0) {
          checkHess <- FALSE
        }
      }
      
      stopflag <- ifelse(checkHess, (fit$output$status[[1]] <= 
          greenOK) & (all(eigen(fit$output$calculatedHessian, 
            symmetric = T, only.values = T)$values > 0)) & 
          (fit$output$minimum <= fit2beat) & (fit$output$minimum <= bestfitsofar), (fit$output$status[[1]] <=  #added bestfitsofar condition
              greenOK) & (fit$output$minimum <= fit2beat) & (fit$output$minimum <= bestfitsofar) )
      if (!stopflag) {
        model <- fit
        newparams<-omxGetParameters(fit) #get recent fit
        if(exists('bestfit')) newparams<-bestfit.params #if bestfit exists use this instead
        if(numdone %% 4 == 0) newparams<-inits #sometimes, use initial start values instead
        #         if(numdone %% 5 == 0) { #sometimes, switch optimizers
        #           if(mxOption(NULL, "Default optimizer")=='CSOLNP') newoptimizer<-'NPSOL'
        #           if(mxOption(NULL, "Default optimizer")=='NPSOL') newoptimizer<-'CSOLNP'
        #           message(paste0('Switching to ',newoptimizer,' optimizer for OpenMx temporarily')) 
        #           mxOption(NULL, "Default optimizer", newoptimizer)
        #         }
        model <- omxSetParameters(model, labels = names(params), 
          values = newparams * runif(length(params),loc-scale,loc+scale))
        fit2beat <- ifelse(fit$output$minimum < fit2beat, fit$output$minimum, 
          fit2beat)
        if(iterationSummary==TRUE){
          message(paste0("Attempt ",numdone," fit:  "))
          message(paste(names(params),": ", fit$output$estimate,"\n"))
          message(paste0("-2LL = ", fit$output$Minus2LogLikelihood))
        }
      }
      
      if(stopflag){
        message('\nSolution found\n')
        fit<-bestfit
        if(length(fit$intervals)>0){ #only calculate confidence intervals once the best fit is established
          fit<-omxSetParameters(fit, labels=names(bestfit.params),values=bestfit.params)
          message("Refit using best inits and estimate confidence intervals\n") 
          #           mxOption(NULL, "Default optimizer", "NPSOL")
          cifit<-suppressWarnings(try(mxRun(fit,intervals=TRUE,suppressWarnings=T,silent=T)))
          if(class(cifit) == "try-error" || cifit$output$status$status== -1) {
            message('Confidence interval estimation generated errors\n')
          } else {
            if (length(summary(cifit)$npsolMessage) > 0) message('Warning messages generated from confidence interval refit\n')
            fit<-cifit
          }
          
        }
        if (length(summary(fit)$npsolMessage) > 0) {
          warning(summary(fit)$npsolMessage)
        }
        
        params <- bestfit.params
        
        if(iterationSummary==TRUE){
          message(paste(names(bestfit.params),": ", bestfit$output$estimate,"\n"))
          message(paste0("-2LL = ", bestfit$output$Minus2LogLikelihood))
        }
        
      }
    } #end 'if fit not an error' section
    if (numdone > extraTries & stopflag==FALSE) { #added stopflag==FALSE
      message('\nRetry limit reached')
      stopflag <- TRUE
      if (exists("bestfit")) {
        fit <- bestfit
        params <- bestfit.params
        if(length(fit$intervals)>0){ #calculate intervals for best fit, even though imperfect
          fit<-omxSetParameters(fit, labels=names(bestfit.params),values=bestfit.params)
          message("Refit using best inits and estimate confidence intervals\n") 
          #           mxOption(NULL, "Default optimizer", "NPSOL")
          cifit<-suppressWarnings(try(mxRun(fit,intervals=TRUE,suppressWarnings=T,silent=T)))
          if(class(cifit) == "try-error" || cifit$output$status$status== -1) {
            message('Confidence interval estimation generated errors, returning fit without confidence intervals\n')
          } else {
            fit<-cifit
          }
        }
        if (length(fit$output$status$statusMsg) > 0) { 
          warning(fit$output$status$statusMsg)
        }
        if(fit$output$status$code==6) message('\nUncertain solution found - consider parameter validity, try again, increase extraTries, change inits, change model, or check data!\n')
        if(iterationSummary==TRUE){
          message(paste(names(bestfit.params),": ", bestfit$output$estimate,"\n"))
          message(paste0("-2LL = ", bestfit$output$Minus2LogLikelihood))
        }
      }
      if (!exists("bestfit")) {
        if (length(fit$output$status$statusMsg) > 0) { 
          warning(fit$output$status$statusMsg)
        }
      }
    }
  }
  if(bestInitsOutput){
    message("\nStart values from best fit:")
    if(paste) message(paste(params, sep=",", collapse = ",")) 
    if(!paste)  message(paste(names(params),": ", params,"\n"))
  }
  return(fit)
}

  
