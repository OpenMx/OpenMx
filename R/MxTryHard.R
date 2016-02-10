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
											initialGradientIterations = as.integer(options()$mxOption$'Gradient iterations'),
											initialTolerance=as.numeric(options()$mxOption$'Optimality tolerance'), 
											checkHess = TRUE, fit2beat = Inf, paste = TRUE,
											iterationSummary=FALSE, bestInitsOutput=TRUE, showInits=FALSE, verbose=0, intervals = FALSE,
											finetuneGradient=TRUE, jitterDistrib=c("rnorm","runif","rcauchy"), exhaustive=FALSE
){
	
	#Initialize stuff:
	jitterDistrib <- match.arg(jitterDistrib)
	if (!is.null(model@compute) && (!.hasSlot(model@compute, '.persist') || !model@compute@.persist)) {
		model@compute <- NULL
	}
	defaultComputePlan <- (is.null(model@compute) || is(model@compute, 'MxComputeDefault'))
	relevantOptions <- list(base::options()$mxOption$"Calculate Hessian", base::options()$mxOption$"Standard Errors")
	if("Calculate Hessian" %in%  names(model@options)){relevantOptions[[1]] <- model@options$"Calculate Hessian"}
	if("Standard Errors" %in%  names(model@options)){relevantOptions[[2]] <- model@options$"Standard Errors"}
	#If the options call for SEs and/or Hessian, there is no custom compute plan, and the Hessian will not be checked
	#every fit attempt, then computing SEs and/or Hessian can be put off until the MLE is obtained:
	SElater <- ifelse( (!checkHess && relevantOptions[[2]]=="Yes" && defaultComputePlan), TRUE, FALSE )
	Hesslater <- ifelse( (!checkHess && relevantOptions[[1]]=="Yes" && defaultComputePlan), TRUE, FALSE )
	doIntervals <- ifelse ( (length(model@intervals) && intervals), TRUE, FALSE )
	lastNoError<-TRUE
	generalTolerance <- 1e-5 #used for hessian check and lowest min check
	gradientStepSize <- initialGradientStepSize
	tolerance <- initialTolerance
	gradientIterations <- initialGradientIterations  
	lastBestFitCount<-0 #number of consecutive improvements in fit
	stopflag <- FALSE #should the iterative optimization process stop
	numdone <- 0
	lowestminsofar<-Inf 
	finalfit<- NULL
	inits<-omxGetParameters(model)
	parlbounds <- omxGetParameters(model=model,fetch="lbound")
	parlbounds[is.na(parlbounds)] <- -Inf
	parubounds <- omxGetParameters(model=model,fetch="ubound")
	parubounds[is.na(parubounds)] <- Inf
	
	
	
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
				#values = params * imxJiggle(dsn=jitterDistrib,n=length(params),loc=loc,scale=scale) + 
				#	imxJiggle(dsn=jitterDistrib,n=length(params),loc=0,scale=scale)
				values=imxJiggle(params=params,lbounds=parlbounds,ubounds=parubounds,dsn=jitterDistrib,loc=loc,scale=scale)
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
				if(lastBestFitCount > 2) model <- omxSetParameters(
					model, labels = names(bestfit.params), 
					#values = params * imxJiggle(dsn=jitterDistrib,n=length(params),loc=loc,scale=scale/10) + 
					#	imxJiggle(dsn=jitterDistrib,n=length(params),loc=0,scale=scale/10)
					values=imxJiggle(params=params,lbounds=parlbounds,ubounds=parubounds,dsn=jitterDistrib,loc=loc,
													 scale=scale/10)
				)
			}
			else{
				model <- omxSetParameters(
					model, labels = names(bestfit.params), 
					#values = params * 
					#	imxJiggle(dsn=jitterDistrib,n=length(params),loc=loc,scale=scale/ifelse(finetuneGradient,10,1)) + 
					#	imxJiggle(dsn=jitterDistrib,n=length(params),loc=0,scale=scale/ifelse(finetuneGradient,10,1))
					values=imxJiggle(params=params,lbounds=parlbounds,ubounds=parubounds,dsn=jitterDistrib,loc=loc,
													 scale=scale/ifelse(finetuneGradient,10,1))
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
				if(any(Hesslater,SElater,doIntervals)){
					message("Running final fit, for Hessian and/or standard errors and/or confidence intervals\n")
					finalfit <- bestfit
					if(defaultComputePlan){
						steps <- list()
						if(doIntervals){
							steps <- c(steps,CI=mxComputeConfidenceInterval(
								plan=mxComputeGradientDescent(
									nudgeZeroStarts=FALSE,gradientIterations=gradientIterations, tolerance=tolerance, 
									maxMajorIter=3000),
								constraintType=ifelse(mxOption(NULL, "Default optimizer") == 'NPSOL','none','ineq')))
						}
						if(Hesslater){steps <- c(steps,ND=mxComputeNumericDeriv())}
						if(SElater){
							steps <- c(steps,SE=mxComputeStandardError(),HQ=mxComputeHessianQuality())
						}
						steps <- c(steps,RD=mxComputeReportDeriv())
						finalfit <- OpenMx::mxModel(finalfit,mxComputeSequence(steps=steps))
					}
					finalfit <- suppressWarnings(try(mxRun(finalfit, suppressWarnings = T, silent=T,	intervals=doIntervals)))
					if(class(finalfit) == "try-error" || finalfit$output$status$status== -1) {
						message('Errors during final fit for Hessian/SEs/CIs\n')
					} else {
						if (length(summary(finalfit)$npsolMessage) > 0) message('Warning messages generated from final fit for final fit for Hessian/SEs/CIs\n')
					}
				}
				if (length(summary(bestfit)$npsolMessage) > 0) {
					warning(summary(bestfit)$npsolMessage)
				}
				if(iterationSummary==TRUE){
					message(paste(names(bestfit.params),": ", bestfit$output$estimate,"\n"))
					message(paste0("-2LL = ", bestfit$output$Minus2LogLikelihood))
				}
				bestfit <- THFrankenmodel(finalfit,bestfit,defaultComputePlan,Hesslater,SElater,doIntervals,checkHess)
			}
		} #end 'if fit not an error' section
		
		
		if (numdone > extraTries && stopflag==FALSE) { #added stopflag==FALSE
			message('\nRetry limit reached')
			stopflag <- TRUE
			if (exists("bestfit")) {
				if(any(Hesslater,SElater,doIntervals)){
					message("Computing Hessian and/or standard errors and/or confidence intervals from imperfect solution\n")
					finalfit <- bestfit
					if(defaultComputePlan){
						steps <- list()
						if(doIntervals){
							steps <- c(steps,CI=mxComputeConfidenceInterval(
								plan=mxComputeGradientDescent(
									nudgeZeroStarts=FALSE,gradientIterations=gradientIterations, tolerance=tolerance, 
									maxMajorIter=3000),
								constraintType=ifelse(mxOption(NULL, "Default optimizer") == 'NPSOL','none','ineq')))
						}
						if(Hesslater){steps <- c(steps,ND=mxComputeNumericDeriv())}
						if(SElater){
							steps <- c(steps,SE=mxComputeStandardError(),HQ=mxComputeHessianQuality())
						}
						steps <- c(steps,RD=mxComputeReportDeriv())
						finalfit <- OpenMx::mxModel(bestfit,mxComputeSequence(steps=steps))
					}
					finalfit <- suppressWarnings(try(mxRun(finalfit, suppressWarnings = T, silent=T,	intervals=doIntervals)))
					if(class(finalfit) == "try-error" || finalfit$output$status$status== -1) {
						message('Errors occurred during final fit for Hessian/SEs/CIs; returning best fit as-is\n')
					}
					if (length(bestfit$output$status$statusMsg) > 0) { 
						warning(bestfit$output$status$statusMsg)
					}
					if(bestfit$output$status$code==6) message('\nUncertain solution found - consider parameter validity, try again, increase extraTries, change inits, change model, or check data!\n')
					if(iterationSummary==TRUE){
						message(paste(names(bestfit.params),": ", bestfit$output$estimate,"\n"))
						message(paste0("-2LL = ", bestfit$output$Minus2LogLikelihood))
					}
					bestfit <- THFrankenmodel(finalfit,bestfit,defaultComputePlan,Hesslater,SElater,doIntervals,checkHess)
	}}}} #end while loop
	
	
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
##' Jiggle parameter values, subject to box constraints.
##' For internal use by mxTryHard().
##' This is an internal function exported for those people who know
##' what they are doing.
##' @param dsn Character string naming the family of distribution to be used.
##' @param n Non-negative integer number of random draws to generate.
##' @param loc Numeric vector of location parameters (medians).
##' @param scale Numeric vector of scale parameters (see source code).
##' @aliases
##' imxJiggle
imxJiggle <- function(params, lbounds, ubounds, dsn, loc, scale){
	if( !(dsn %in% c("rnorm","runif","rcauchy")) ){stop("unrecognized value for argument 'dsn'")}
	n <- length(params)
	if(dsn=="rnorm"){
		out <- params * rnorm(n=n,mean=loc,sd=scale)) + rnorm(n=n,mean=0,sd=scale)
	}
	if(dsn=="runif"){
		out <- params * runif(n=n,min=loc-scale,max=loc+scale) + runif(n=n,min=0-scale,max=scale)
	}
	if(dsn=="rcauchy"){
		out <- params * rcauchy(n=n,location=loc,scale=scale) + rcauchy(n=n,location=0,scale=scale)
	}
	if(any(out<lbounds)){out[out<lbounds] <- lbounds[out<lbounds]}
	if(any(out>ubounds)){out[out>ubounds] <- ubounds[out>ubounds]}
	return(out)
}



THFrankenmodel <- function(finalfit,bestfit,defaultComputePlan,Hesslater,SElater,doIntervals,checkHess){
	if( is.null(finalfit) || !any(Hesslater,SElater,doIntervals) || ("try-error" %in% class(finalfit)) || 
			finalfit$output$status$status== -1 ){return(bestfit)}
	if(defaultComputePlan){
		steps <- list(GD=bestfit@compute@steps[["GD"]])
		if(doIntervals){steps <- c(steps,CI=finalfit@compute@steps[["CI"]])}
		if(Hesslater || SElater){
			if(Hesslater){steps <- c(steps,ND=finalfit@compute@steps[["ND"]])}
			if(SElater){steps <- 	c(steps,SE=finalfit@compute@steps[["SE"]],HQ=finalfit@compute@steps[["HQ"]])}
			steps <- c(steps,RD=finalfit@compute@steps[["RD"]])
		}
		else{
			if(checkHess){steps <- c(steps,ND=bestfit@compute@steps[["ND"]],SE=bestfit@compute@steps[["SE"]])}
			steps <- c(steps,RD=bestfit@compute@steps[["RD"]])
		}
		steps <- c(steps,RE=bestfit@compute@steps[["RE"]])
		bestfit@compute@steps <- steps
	}
	else{
		if( doIntervals && ("MxComputeConfidenceInterval" %in% unlist(lapply(bestfit@compute@steps,class))) &&
			 ("MxComputeConfidenceInterval" %in% unlist(lapply(finalfit@compute@steps,class))) ){
			f <- which("MxComputeConfidenceInterval"==unlist(lapply(finalfit@compute@steps,class)))
			t <- which("MxComputeConfidenceInterval"==unlist(lapply(bestfit@compute@steps,class)))
			bestfit@compute@steps[t] <- finalfit@compute@steps[f]
	}}
	bestfit@output$timestamp <- finalfit@output$timestamp
	if(doIntervals){
		bestfit@output$confidenceIntervals <- finalfit@output$confidenceIntervals
		bestfit@output$confidenceIntervalCodes <- finalfit@output$confidenceIntervalCodes
	}
	if(Hesslater || SElater){
		bestfit@output$calculatedHessian <- finalfit@output$calculatedHessian
		bestfit@output$hessian <- finalfit@output$hessian
		bestfit@output$standardErrors <- finalfit@output$standardErrors
		bestfit@output$infoDefinite <- finalfit@output$infoDefinite
		bestfit@output$conditionNumber <- finalfit@output$conditionNumber
	}
	bestfit@output$evaluations <- bestfit@output$evaluations + finalfit@output$evaluations
	bestfit@output$frontendTime <- bestfit@output$frontendTime + finalfit@output$frontendTime
	bestfit@output$backendTime <- bestfit@output$backendTime + finalfit@output$backendTime
	bestfit@output$independentTime <- bestfit@output$independentTime + finalfit@output$independentTime
	bestfit@output$wallTime <- bestfit@output$wallTime + finalfit@output$wallTime
	bestfit@output$cpuTime <- bestfit@output$cpuTime + finalfit@output$cpuTime
	bestfit@.modifiedSinceRun <- FALSE
	return(bestfit)
}

