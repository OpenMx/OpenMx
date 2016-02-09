#
#   Copyright 2007-2016 The OpenMx Project
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

mxTryHard <- function(model, extraTries = 10, greenOK = FALSE, loc = 1, 
											scale = 0.25,  initialGradientStepSize = .00001, 
											initialGradientIterations = options()$mxOption$'Gradient iterations',
											initialTolerance=options()$mxOption$'Optimality tolerance', 
											checkHess = TRUE, fit2beat = Inf, paste = TRUE,
											iterationSummary=FALSE, bestInitsOutput=TRUE, showInits=FALSE, verbose=0, intervals = FALSE,
											finetuneGradient=TRUE, jitterDistrib=c("rnorm","runif","rcauchy"), exhaustive=FALSE
){
	
	#Initialize stuff:
	jitterDistrib <- match.call()$jitterDistrib[1]
	relevantOptions <- list(base::options()$mxOption$"Calculate Hessian", base::options()$mxOption$"Standard Errors")
	if("Calculate Hessian" %in%  names(model@options)){relevantOptions[[1]] <- model@options$"Calculate Hessian"}
	if("Standard Errors" %in%  names(model@options)){relevantOptions[[2]] <- model@options$"Standard Errors"}
	if (!is.null(model@compute) && (!.hasSlot(model@compute, '.persist') || !model@compute@.persist)) {
		model@compute <- NULL
	}
	defaultComputePlan <- (is.null(model@compute) || is(model@compute, 'MxComputeDefault'))
	lastNoError<-TRUE
	generalTolerance <- 1e-5 #used for hessian check and lowest min check
	gradientStepSize <- initialGradientStepSize
	tolerance <- initialTolerance
	gradientIterations <- initialGradientIterations  
	lastBestFitCount<-0 #number of consecutive improvements in fit
	stopflag <- FALSE #should the iterative optimization process stop
	numdone <- 0
	lowestminsofar<-Inf 
	cifit <- sefit <- NULL
	inits<-omxGetParameters(model) 
	
	
	#Begin main loop.
	while (!stopflag) {
		message(paste0('\nBegin fit attempt ', numdone+1, ' of at maximum ', extraTries +1, ' tries'))
		if(lastNoError==TRUE) params <- omxGetParameters(model)
		
		
		if(lastBestFitCount == 0 && numdone > 0){ #if the last fit was not the best
			if(exists('bestfit')) params <- bestfit.params #if bestfit exists use this instead
			if(numdone %% 4 == 0 && finetuneGradient) params <- inits #sometimes, use initial start values instead
			#^^^No reason to re-use start values unless optimization-control parameters have been changed,
			#which will only be happening when finetuneGradient is TRUE.
			model <- omxSetParameters(
				model, labels = names(params), 
				values = params * imxJiggle(dsn=jitterDistrib,n=length(params),loc=loc,scale=scale) + 
					imxJiggle(dsn=jitterDistrib,n=length(params),loc=0,scale=scale)
			)
			if(finetuneGradient){
				gradientStepSize <- initialGradientStepSize
				tolerance <- initialTolerance
				gradientIterations<-initialGradientIterations
			}
		}#end if last fit not best section
		
		
		if(lastBestFitCount > 0){ #if the last fit was the best so far
			if(exists('bestfit')) {
				params <- bestfit.params      
				model <- bestfit
			}
			if(defaultComputePlan==TRUE && finetuneGradient){
				if(lastBestFitCount == 2) gradientStepSize <- gradientStepSize *.1
				if(lastBestFitCount == 3) gradientStepSize <- gradientStepSize *10
				if(lastBestFitCount == 5) gradientStepSize <- gradientStepSize *.1
				if(lastBestFitCount  > 0) tolerance<-tolerance * .001 
				if(lastBestFitCount  > 0) gradientIterations<-gradientIterations+2
				# if(lastBestFitCount  %in% seq(4,100,4)) gradientIterations<-gradientIterations-1
				if(lastBestFitCount > 2) model <- omxSetParameters(
					model, labels = names(bestfit.params), 
					values = params * imxJiggle(dsn=jitterDistrib,n=length(params),loc=loc,scale=scale/10) + 
						imxJiggle(dsn=jitterDistrib,n=length(params),loc=0,scale=scale/10)
				)
			}
			else{
				model <- omxSetParameters(
					model, labels = names(bestfit.params), 
					values = params * 
						imxJiggle(dsn=jitterDistrib,n=length(params),loc=loc,scale=scale/ifelse(finetuneGradient,10,1)) + 
						imxJiggle(dsn=jitterDistrib,n=length(params),loc=0,scale=scale/ifelse(finetuneGradient,10,1))
				)
			}
		}#end if last fit was best section
		
		
		if(defaultComputePlan==TRUE){
			steps <- list(GD=mxComputeGradientDescent(
				verbose=verbose, gradientStepSize = gradientStepSize, 
				nudgeZeroStarts=FALSE,   gradientIterations = gradientIterations, tolerance=tolerance, 
				maxMajorIter=3000))
			if(checkHess){steps <- c(steps,ND=mxComputeNumericDeriv(),SE=mxComputeStandardError())}
			model <- OpenMx::mxModel(
				model,
				mxComputeSequence(c( steps,RD=mxComputeReportDeriv(),RE=mxComputeReportExpectation() )))
			rm(steps)
		}
		if(showInits==TRUE) {
			message('Starting values:  ')
			message(paste0(names(omxGetParameters(model)),' : ', omxGetParameters(model),'\n'))
		}
		fit <- suppressWarnings(try(mxRun(model, suppressWarnings = T, unsafe=T, silent=T,intervals=FALSE)))
		numdone <- numdone + 1
		if( class(fit) == "try-error" || is.na(fit$output$minimum) || fit$output$status$status== -1) {
			lastBestFitCount <- 0
			lastNoError<-FALSE
			message('\n Fit attempt generated errors') 
		}
		
		
		if(class(fit) != "try-error" && !is.na(fit$output$minimum) && fit$output$status$status != -1){ #if fit was not an error
			if (fit$output$minimum >= lowestminsofar + generalTolerance) {
				lastBestFitCount <- 0
				lastNoError<-TRUE
				message(paste0('\n Fit attempt worse than current best:  ',fit$output$minimum ,' vs ', lowestminsofar )) 
			}
			if(fit$output$minimum >= lowestminsofar) lastBestFitCount<-0
			if (fit$output$minimum < lowestminsofar && is.finite(fit$output$minimum)) { #if this is the best fit so far
				message(paste0('\n Lowest minimum so far:  ',fit$output$minimum) )
				lastBestFitCount<-lastBestFitCount+1 
				lowestminsofar <- fit$output$minimum
				lastNoError<-TRUE
				bestfit <- fit
				bestfit.params <- omxGetParameters(bestfit)
			}
			if (fit$output$minimum <= lowestminsofar + generalTolerance) { #if this is the best fit or equal best so far, check the following
				###########stopflag checks
				stopflag<-TRUE
				if(fit$output$status[[1]] > greenOK) {
					stopflag<-FALSE
					message(paste0('\n OpenMx status code ', fit$output$status[[1]], ' greater than ', as.numeric(greenOK)))
				}
				if(fit$output$minimum > fit2beat) {
					message(paste0('\n Fit value of ', fit$output$minimum, ' greater than fit2beat of ', fit2beat))
					stopflag<-FALSE
				}
				if(fit$output$minimum > lowestminsofar + generalTolerance){
					message(paste0('\n Fit value of ', fit$output$minimum, ' greater than best so far of ', lowestminsofar))
					stopflag<-FALSE
				}
				if(checkHess==TRUE) {
					hessEigenval <- try(eigen(fit$output$calculatedHessian, symmetric = T, only.values = T)$values)
					if(class(hessEigenval)=='try-error') {
						message(paste0('\n Eigenvalues of hessian could not be calculated'))
						stopflag<-FALSE
					}
					if(class(hessEigenval)!='try-error' && any(hessEigenval < 0)) {
						message(paste0('\n Not all eigenvalues of hessian are greater than ', 0,': ', paste(hessEigenval,collapse=', ')))
						stopflag<-FALSE
					}
					if(stopflag ==TRUE) bestfit <- fit
				}
			}#end stopflag checks and if lowest min section
			
			if (!stopflag){
				if(iterationSummary==TRUE){
					message(paste0("\n Attempt ",numdone," fit:  "))
					message(paste(names(params),": ", fit$output$estimate,"\n"))
					message(paste0("-2LL = ", fit$output$Minus2LogLikelihood))
				}
			}
			
			if(stopflag){
				message('\nSolution found\n')
				if( !checkHess && 
						(relevantOptions$"Calculate Hessian"=="Yes" || relevantOptions$"Standard Errors"=="Yes") ){
					message("Computing standard errors and/or explicitly-calculated Hessian\n")
					steps <- list()
					if(relevantOptions$"Calculate Hessian"=="Yes"){steps <- c(steps,ND=mxComputeNumericDeriv())}
					if(relevantOptions$"Standard Errors"=="Yes"){
						steps <- c(steps,SE=mxComputeStandardError(),HQ=mxComputeHessianQuality())
					}
					steps <- c(steps,RD=mxComputeReportDeriv())
					bestfit <- OpenMx::mxModel(bestfit,mxComputeSequence(steps=steps))
					sefit <- suppressWarnings(try(mxRun(bestfit, suppressWarnings = T, unsafe=T, silent=T,intervals=FALSE)))
					rm(steps)
					if(class(sefit) == "try-error" || sefit$output$status$status== -1) {
						message('Errors during SE/Hessian computation\n')
					} else {
						if (length(summary(sefit)$npsolMessage) > 0) message('Warning messages generated from SE/Hessian refit\n')
					}
				}
				if(length(bestfit$intervals)>0 && intervals==TRUE){ #only calculate confidence intervals once the best fit is established
					message("Estimating confidence intervals\n") 
					if(defaultComputePlan==TRUE) bestfit <- OpenMx::mxModel(
						bestfit, 
						mxComputeSequence(list(
							CI=mxComputeConfidenceInterval(
								plan=mxComputeGradientDescent(
									nudgeZeroStarts=FALSE,gradientIterations=gradientIterations, tolerance=tolerance, 
									maxMajorIter=3000),
								constraintType=ifelse(mxOption(NULL, "Default optimizer") == 'NPSOL','none','ineq'))
							#mxComputeNumericDeriv(), mxComputeStandardError(), 
							#mxComputeReportDeriv())))
						)))
					cifit<-suppressWarnings(try(mxRun(bestfit,intervals=TRUE,suppressWarnings=T,silent=T)))
					if(class(cifit) == "try-error" || cifit$output$status$status== -1) {
						message('Confidence interval estimation generated errors\n')
					} else {
						if (length(summary(cifit)$npsolMessage) > 0) message('Warning messages generated from confidence interval refit\n')
					}
				}
				if (length(summary(bestfit)$npsolMessage) > 0) {
					warning(summary(bestfit)$npsolMessage)
				}
				if(iterationSummary==TRUE){
					message(paste(names(bestfit.params),": ", bestfit$output$estimate,"\n"))
					message(paste0("-2LL = ", bestfit$output$Minus2LogLikelihood))
				}
				bestfit <- THFrankenmodel(cifit,sefit,bestfit,defaultComputePlan)
			}
		} #end 'if fit not an error' section
		
		
		if (numdone > extraTries && stopflag==FALSE) { #added stopflag==FALSE
			message('\nRetry limit reached')
			stopflag <- TRUE
			if (exists("bestfit")) {
				if( !checkHess && 
						(relevantOptions$"Calculate Hessian"=="Yes" || relevantOptions$"Standard Errors"=="Yes") ){
					message("Computing SEs and/or Hessian for imperfect solution\n")
					steps <- list()
					if(relevantOptions$"Calculate Hessian"=="Yes"){steps <- c(steps,ND=mxComputeNumericDeriv())}
					if(relevantOptions$"Standard Errors"=="Yes"){
						steps <- c(steps,SE=mxComputeStandardError(),HQ=mxComputeHessianQuality())
					}
					steps <- c(steps,RD=mxComputeReportDeriv())
					bestfit <- OpenMx::mxModel(bestfit,mxComputeSequence(steps=steps))
					sefit <- suppressWarnings(try(mxRun(bestfit, suppressWarnings = T, unsafe=T, silent=T,intervals=FALSE)))
					rm(steps)
					if(class(sefit) == "try-error" || sefit$output$status$status== -1) {
						message('Errors during SE/Hessian computation; returning fit without them\n')
					}
				}
				if(length(bestfit$intervals)>0 && intervals==TRUE){ #calculate intervals for best fit, even though imperfect
					message("Estimating confidence intervals for imperfect solution\n") 
					if(defaultComputePlan==TRUE) bestfit <- OpenMx::mxModel(
						bestfit, 
						mxComputeSequence(list(
							mxComputeConfidenceInterval(
								plan=mxComputeGradientDescent(
									nudgeZeroStarts=FALSE, 
									gradientIterations=gradientIterations, tolerance=tolerance, 
									maxMajorIter=3000),
								constraintType=ifelse(mxOption(NULL, "Default optimizer") == 'NPSOL','none','ineq'))
							#mxComputeNumericDeriv(), mxComputeStandardError(), 
							#mxComputeReportDeriv())))
						)))
					
					cifit<-suppressWarnings(try(mxRun(bestfit,intervals=TRUE,suppressWarnings=T,silent=T)))
					if(class(cifit) == "try-error" || cifit$output$status$status== -1) {
						message('Confidence interval estimation generated errors, returning fit without confidence intervals\n')
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
				bestfit <- THFrankenmodel(cifit,sefit,bestfit,defaultComputePlan)
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
	
	if( defaultComputePlan && !("try-error" %in% class(bestfit)) ){bestfit@compute@.persist <- FALSE}
	
	return(bestfit)
}



##' imxJiggle
##' 
##' Wrapper to random-number generators.
##' This is an internal function exported for those people who know
##' what they are doing.
##' @param dsn Character string naming the family of distribution to be used.
##' @param n Non-negative integer number of random draws to generate.
##' @param loc Numeric vector of location parameters (medians).
##' @param scale Numeric vector of scale parameters (see source code).
##' @aliases
##' imxJiggle
imxrunif <- function(dsn, n, loc, scale){
	if( !(dsn %in% c("rnorm","runif","rcauchy")) ){stop("unrecognized value for argument 'dsn'")}
	if(dsn=="rnorm"){return(rnorm(n=n,mean=loc,sd=scale))}
	if(dsn=="runif"){return(runif(n=n,min=loc-scale,max=loc+scale))}
	if(dsn=="rcauchy"){return(rcauchy(n=n,location=loc,scale=scale))}
}



THFrankenmodel <- function(cifit,sefit,bestfit,defaultComputePlan){
	if(is.null(cifit) && is.null(sefit)){return(bestfit)}
	if(!is.null(cifit) && is.null(sefit)){
		if(defaultComputePlan){
			bestfit@compute@steps <- list(
				GD=bestfit@compute@steps[["GD"]],ND=bestfit@compute@steps[["ND"]],
				SE=bestfit@compute@steps[["SE"]],RD=bestfit@compute@steps[["RD"]],
				RE=bestfit@compute@steps[["RE"]],CI=cifit@compute@steps[[1]])
		}
		else{bestfit@compute@steps[[length(bestfit@compute@steps)+1]] <- cifit@compute@steps[[1]]}
		bestfit@output$confidenceIntervals <- cifit@output$confidenceIntervals
		bestfit@output$confidenceIntervalCodes <- cifit@output$confidenceIntervalCodes
		bestfit@output$timestamp <- cifit@output$timestamp
		bestfit@output$evaluations <- bestfit@output$evaluations + cifit@output$evaluations
		bestfit@output$frontendTime <- bestfit@output$frontendTime + cifit@output$frontendTime
		bestfit@output$backendTime <- bestfit@output$backendTime + cifit@output$backendTime
		bestfit@output$independentTime <- bestfit@output$independentTime + cifit@output$independentTime
		bestfit@output$wallTime <- bestfit@output$wallTime + cifit@output$wallTime
		bestfit@output$cpuTime <- bestfit@output$cpuTime + cifit@output$cpuTime
		bestfit@.modifiedSinceRun <- FALSE
		return(bestfit)
	}
	if(is.null(cifit) && !is.null(sefit)){
		if(defaultComputePlan){
			bestfit@compute@steps <- list(
				GD=bestfit@compute@steps[["GD"]],ND=sefit@compute@steps[["ND"]],
				SE=sefit@compute@steps[["SE"]],HQ=sefit@compute@steps[["HQ"]],
				RD=sefit@compute@steps[["RD"]],RE=bestfit@compute@steps[["RE"]])
		}
		else{bestfit@compute@steps <- c(bestfit@compute@steps, sefit@compute@steps)}
		bestfit@output$calculatedHessian <- sefit@output$calculatedHessian
		bestfit@output$hessian <- sefit@output$hessian
		bestfit@output$standardErrors <- sefit@output$standardErrors
		bestfit@output$infoDefinite <- sefit@output$infoDefinite
		bestfit@output$conditionNumber <- sefit@output$conditionNumber
		bestfit@output$timestamp <- sefit@output$timestamp
		bestfit@output$evaluations <- bestfit@output$evaluations + sefit@output$evaluations
		bestfit@output$frontendTime <- bestfit@output$frontendTime + sefit@output$frontendTime
		bestfit@output$backendTime <- bestfit@output$backendTime + sefit@output$backendTime
		bestfit@output$independentTime <- bestfit@output$independentTime + sefit@output$independentTime
		bestfit@output$wallTime <- bestfit@output$wallTime + sefit@output$wallTime
		bestfit@output$cpuTime <- bestfit@output$cpuTime + sefit@output$cpuTime
		bestfit@.modifiedSinceRun <- FALSE
		return(bestfit)
	}
	if(!is.null(cifit) && !is.null(sefit)){
		if(defaultComputePlan){
		bestfit@compute@steps <- list(
			GD=bestfit@compute@steps[["GD"]],ND=sefit@compute@steps[["ND"]],
			SE=sefit@compute@steps[["SE"]],HQ=sefit@compute@steps[["HQ"]],
			RD=sefit@compute@steps[["RD"]],RE=bestfit@compute@steps[["RE"]],CI=cifit@compute@steps[[1]])
		}
		else{bestfit@compute@steps <- c(bestfit@compute@steps, sefit@compute@steps, cifit@compute@steps}
		bestfit@output$confidenceIntervals <- cifit@output$confidenceIntervals
		bestfit@output$confidenceIntervalCodes <- cifit@output$confidenceIntervalCodes
		bestfit@output$calculatedHessian <- sefit@output$calculatedHessian
		bestfit@output$hessian <- sefit@output$hessian
		bestfit@output$standardErrors <- sefit@output$standardErrors
		bestfit@output$infoDefinite <- sefit@output$infoDefinite
		bestfit@output$conditionNumber <- sefit@output$conditionNumber
		bestfit@output$timestamp <- cifit@output$timestamp
		bestfit@output$evaluations <- bestfit@output$evaluations + sefit@output$evaluations + cifit@output$evaluations
		bestfit@output$frontendTime <- bestfit@output$frontendTime + sefit@output$frontendTime + 
			cifit@output$frontendTime
		bestfit@output$backendTime <- bestfit@output$backendTime + sefit@output$backendTime + cifit@output$backendTime
		bestfit@output$independentTime <- bestfit@output$independentTime + sefit@output$independentTime + 
			cifit@output$independentTime
		bestfit@output$wallTime <- bestfit@output$wallTime + sefit@output$wallTime + cifit@output$wallTime
		bestfit@output$cpuTime <- bestfit@output$cpuTime + sefit@output$cpuTime + cifit@output$cpuTime
		bestfit@.modifiedSinceRun <- FALSE
		return(bestfit)
	}
}

