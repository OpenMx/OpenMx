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

mxRun <- function(model, ..., intervals=NULL, silent = FALSE, 
		suppressWarnings = FALSE, unsafe = FALSE,
		checkpoint = FALSE, useSocket = FALSE, onlyFrontend = FALSE, 
		useOptimizer = TRUE){

	if (is.null(intervals)) {
		# OK
	} else if (length(intervals) != 1 ||
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
		compute <- NULL
		fitNum <- paste(model@name, 'fitfunction', sep=".")
		if (!useOptimizer) {
			compute <- mxComputeOnce(from=fitNum, 'fit', .is.bestfit=TRUE)
		} else {
			steps = list(GD=mxComputeGradientDescent(fitfunction=fitNum))
			if (length(intervals) && intervals) {
				ciOpt <- mxComputeGradientDescent(
				    fitfunction=fitNum, nudgeZeroStarts=FALSE, maxMajorIter=150)
				cType <- 'ineq'
				if (ciOpt$engine == "NPSOL") cType <- 'none'
				steps <- c(steps, CI=mxComputeConfidenceInterval(
				    fitfunction=fitNum, constraintType=cType, plan=ciOpt))
			}
			if (options[["Calculate Hessian"]] == "Yes") {
				steps <- c(steps, ND=mxComputeNumericDeriv(fitfunction=fitNum))
			}
			if (options[["Standard Errors"]] == "Yes") {
				steps <- c(steps, SE=mxComputeStandardError(), HQ=mxComputeHessianQuality())
			}
			compute <- mxComputeSequence(c(steps, RD=mxComputeReportDeriv()))
		}
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
		prec <- lapply(expectations, genericExpGetPrecision)

		functionPrecision <- Reduce(max, c(as.numeric(options[['Function precision']]),
						   sapply(prec, function(x) x[['functionPrecision']])))
		options[['Function precision']] <- as.character(functionPrecision)

		if (defaultComputePlan && is(model@compute, "MxComputeSequence")) {
			iterations <- Reduce(min, c(4L, sapply(prec, function(x) x[['iterations']])))
			stepSize <- Reduce(max, c(1e-4, sapply(prec, function(x) x[['stepSize']])))
			model <- adjustDefaultNumericDeriv(model, iterations, stepSize)
			flatModel <- adjustDefaultNumericDeriv(flatModel, iterations, stepSize)
		}
	}

	fitfunctions <- convertFitFunctions(flatModel, model, labelsData, dependencies)
	data <- convertDatasets(flatModel@datasets, model, flatModel)
	numAlgebras <- length(algebras)
	algebras <- append(algebras, fitfunctions)
	constraints <- convertConstraints(flatModel)
	parameters <- flatModel@parameters
	numParam <- length(parameters)
	if (numParam == 0 && defaultComputePlan && !is.null(model@fitfunction)) {
		compute <- mxComputeOnce(from=paste(model@name, 'fitfunction', sep="."),
					 'fit', .is.bestfit=TRUE)
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
	if (onlyFrontend) return(model)
	output <- .Call(backend,
			constraints, matrices, parameters,
			algebras, expectations, computes,
			data, intervalList, communication, options, defVars, PACKAGE = "OpenMx")
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
		names(parameters), output)
	if(length(model$data) > 0 && model$data$type=="acov"){
		wlsSEs <- imxWlsStandardErrors(model)
		model@output$standardErrors <- wlsSEs$SE
		model@output$hessian <- 2*solve(wlsSEs$Cov) #puts in same units as m2ll Hessian
		model@output$calculatedHessian <- model@output$hessian
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
    scale = 0.25,  initialGradientStepSize = .00001, initialGradientIterations = 1,
    initialTolerance=1e-12, checkHess = TRUE, fit2beat = Inf, paste = TRUE,
    iterationSummary=FALSE, bestInitsOutput=TRUE, showInits=FALSE, verbose=0,
    ...){
    lastNoError<-TRUE
    
    gradientStepSize<- initialGradientStepSize
    tolerance <- initialTolerance
    lastBestFitCount<-0 #number of consecutive improvements in fit
    gradientIterations <- initialGradientIterations  
    stopflag <- FALSE #should the iterative optimization process stop
    numdone <- 0
    lowestminsofar<-Inf 
    inits<-omxGetParameters(model) 
    
    
    while (!stopflag) {
      
      if(iterationSummary==TRUE) message(paste0('\nBegin fit attempt ', numdone+1, ' of at maximum ', extraTries +1, ' tries'))
      
      if(lastNoError==TRUE) params <- omxGetParameters(model)
      
      if(lastBestFitCount == 0 & numdone > 0){ #if the last fit was not the best
        message('\n Setting gradientStepSize and gradientIterations to defaults')
        if(exists('bestfit')) params<-bestfit.params #if bestfit exists use this instead
        if(numdone %% 4 == 0) params<-inits #sometimes, use initial start values instead
        #       if(numdone %% 5 == 0) { #sometimes, switch optimizers
        #         if(mxOption(NULL, "Default optimizer")=='CSOLNP') newoptimizer<-'NPSOL'
        #         if(mxOption(NULL, "Default optimizer")=='NPSOL') newoptimizer<-'CSOLNP'
        #         message(paste0('Switching to ',newoptimizer,' optimizer for OpenMx temporarily')) 
        #         mxOption(NULL, "Default optimizer", newoptimizer)
        #       }
        
        model <- omxSetParameters(model, labels = names(params), 
          values = params * runif(length(params),loc-(rnorm(1,scale,1)),loc+(rnorm(1,scale,1))))
        
        gradientStepSize <- initialGradientStepSize
        #             tolerance <- initialTolerance
        gradientIterations<-initialGradientIterations
      }
      
      
      
      
      if(lastBestFitCount > 0){ #if the last fit was the best so far
        if(exists('bestfit')) {
          params<-bestfit.params      
          model<-bestfit
        }
        message('\n Lowest minimum so far')
        
        
        if(lastBestFitCount == 3) gradientStepSize <- gradientStepSize *.1
        if(lastBestFitCount == 4) gradientStepSize <- gradientStepSize *10
        if(lastBestFitCount  %in% seq(1,100,4)) tolerance<-tolerance * .01 
        if(lastBestFitCount  %in% seq(1,100,4)) gradientIterations<-gradientIterations+1
        if(lastBestFitCount  %in% seq(3,100,4)) gradientIterations<-gradientIterations-1
        if(lastBestFitCount > 5) model <- omxSetParameters(model, labels = names(bestfit.params), 
          values = bestfit.params * runif(length(bestfit.params),loc-(rnorm(1,scale,1)/10),loc+(rnorm(1,scale,1)/10)))
      }
      
      
      
      
      model<-mxModel(model, mxComputeSequence(list(
        mxComputeGradientDescent(verbose=verbose, gradientStepSize = gradientStepSize, 
          nudgeZeroStarts=FALSE,   gradientIterations = gradientIterations, tolerance=tolerance, 
          maxMajorIter=3000),
        mxComputeNumericDeriv(), mxComputeStandardError(),  
        mxComputeReportDeriv())))
      
      if(showInits==TRUE) {
        message('Starting values:  ')
        message(paste0(names(omxGetParameters(model)),' : ', omxGetParameters(model),'\n'))
      }
      
      
      fit <- suppressWarnings(try(mxRun(model, suppressWarnings = T, unsafe=T, silent=T,intervals=FALSE)))
      
      numdone <- numdone + 1
      
      
      if(class(fit) == "try-error" || fit$output$status$status== -1 || is.na(fit$output$minimum)) {
        lastBestFitCount <- 0
        lastNoError<-FALSE
        message('\n Fit attempt generated errors') 
      }
      
      if(class(fit) != "try-error" && fit$output$status$status != -1 && !is.na(fit$output$minimum)) { #if fit was not an error
        
        if (fit$output$minimum >= lowestminsofar) {
          lastBestFitCount <- 0
          lastNoError<-TRUE
          message('\n Fit attempt did not beat the current best') 
        }
        
        if (fit$output$minimum < lowestminsofar) { #if this is the best fit so far, check the following if required
          lastBestFitCount<-lastBestFitCount+1 
          lowestminsofar <- fit$output$minimum
          lastNoError<-TRUE
          bestfit <- fit
          bestfit.params <- omxGetParameters(bestfit)
          
          if (length(fit$output$calculatedHessian) == 0) {
            checkHess <- FALSE
          }
          if (checkHess) {
            if (sum(is.na(fit$output$calculatedHessian)) > 
                0) {
              checkHess <- FALSE
            }
          }
        }
        
        
        stopflag <- ifelse(checkHess, (fit$output$status[[1]] <= 
            greenOK) & (all(eigen(fit$output$calculatedHessian, 
              symmetric = T, only.values = T)$values > 0)) & 
            (fit$output$minimum <= fit2beat) & (fit$output$minimum <= lowestminsofar), (fit$output$status[[1]] <=  #added lowestminsofar condition
                greenOK) & (fit$output$minimum <= fit2beat) & (fit$output$minimum <= lowestminsofar) )
        
        if (!stopflag) {        
          if(iterationSummary==TRUE){
            message(paste0("Attempt ",numdone," fit:  "))
            message(paste(names(params),": ", fit$output$estimate,"\n"))
            message(paste0("-2LL = ", fit$output$Minus2LogLikelihood))
          }
        }
        
        if(stopflag){
          message('\nSolution found\n')
          if(length(bestfit$intervals)>0){ #only calculate confidence intervals once the best fit is established
            
            message("Estimating confidence intervals\n") 
            
            mxOption(NULL, "Default optimizer")
            bestfit <- OpenMx::mxModel(bestfit, 
              mxComputeSequence(list(
                mxComputeConfidenceInterval(plan=mxComputeGradientDescent(nudgeZeroStarts=FALSE, 
                  gradientIterations=gradientIterations, tolerance=tolerance, 
                  maxMajorIter=1000),
                  constraintType=ifelse(mxOption(NULL, "Default optimizer") == 'NPSOL','none','ineq')),
                mxComputeNumericDeriv(), mxComputeStandardError(), 
                mxComputeReportDeriv())))
            
            cifit<-suppressWarnings(try(mxRun(bestfit,intervals=TRUE,suppressWarnings=T,silent=T)))
            
            
            if(class(cifit) == "try-error" || cifit$output$status$status== -1) {
              message('Confidence interval estimation generated errors\n')
            } else {
              if (length(OpenMx::summary(cifit)$npsolMessage) > 0) message('Warning messages generated from confidence interval refit\n')
              bestfit<-cifit
            }
            
          }
          if (length(OpenMx::summary(bestfit)$npsolMessage) > 0) {
            warning(OpenMx::summary(bestfit)$npsolMessage)
          }
          
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
          
          if(length(bestfit$intervals)>0){ #calculate intervals for best fit, even though imperfect
            message("Estimate confidence intervals for imperfect solution\n") 
            #           mxOption(NULL, "Default optimizer", "NPSOL")
            cifit<-suppressWarnings(try(mxRun(bestfit,intervals=TRUE,suppressWarnings=T,silent=T)))
            if(class(cifit) == "try-error" || cifit$output$status$status== -1) {
              message('Confidence interval estimation generated errors, returning fit without confidence intervals\n')
            } else {
              bestfit<-cifit
            }
          }
          if (length(bestfit$output$status$statusMsg) > 0) { 
            warning(bestfit$output$status$statusMsg)
          }
          if(bestfit$output$status$code==6) message('\nUncertain solution found - consider parameter validity, try again, increase extraTries, change inits, change model, or check data!\n')
          if(iterationSummary==TRUE){
            message(paste(names(bestfit.params),": ", bestfit$output$estimate,"\n"))
            message(paste0("-2LL = ", bestfit$output$Minus2LogLikelihood))
          }
        }
      }
    } #end while loop
    
    
    
    if(bestInitsOutput && exists("bestfit")){
      bestfit.params <- omxGetParameters(bestfit)
      message("\nStart values from best fit:")
      if(paste) message(paste(bestfit.params, sep=",", collapse = ",")) 
      if(!paste)  message(paste(names(bestfit.params),": ", bestfit.params,"\n"))
    }
    
    if (!exists("bestfit")) {
      if(class(fit) == 'try-error') warning(fit[[length(fit)]])
      message('All fit attempts resulted in errors - check starting values or model specification')
      bestfit<-fit
    }
    
    return(bestfit)
  }
