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

#TODO:
#Need more input checking?  For instance, initialGradientIterations should be a positive integer, right?
#How can mxTryHard() be improved for ordinal-threshold analyses?

mxTryHard <- function(model, extraTries = 10, greenOK = FALSE, loc = 1, 
											scale = 0.25,  initialGradientStepSize = .00001, 
											initialGradientIterations = as.integer(options()$mxOption$'Gradient iterations'),
											initialTolerance=as.numeric(options()$mxOption$'Optimality tolerance'), 
											checkHess = TRUE, fit2beat = Inf, paste = TRUE,
											iterationSummary=FALSE, bestInitsOutput=TRUE, showInits=FALSE, verbose=0, intervals = FALSE,
											finetuneGradient=TRUE, jitterDistrib=c("runif","rnorm","rcauchy"), exhaustive=FALSE,
											maxMajorIter=3000, OKstatuscodes, wtgcsv=c("prev","best","initial")
){
	
	#Initialize stuff & check inputs:
	jitterDistrib <- match.arg(jitterDistrib)
	wtgcsv <- match.arg(wtgcsv,c("prev","best","initial"),several.ok=T)
	if(missing(OKstatuscodes)){OKstatuscodes <- as.integer(c(0,as.logical(greenOK[1])))}
	else if( !(0 %in% OKstatuscodes) ){OKstatuscodes <- c(OKstatuscodes,0)}
	#if( !("MxModel" %in% class(model)) ){stop("argument 'model' must be an object of class 'MxModel'")}
	if(initialTolerance<0){stop("value for argument 'initialTolerance' cannot be negative")}
	if (!is.null(model@compute) && (!.hasSlot(model@compute, '.persist') || !model@compute@.persist)) {
		model@compute <- NULL
	}
	defaultComputePlan <- (is.null(model@compute) || is(model@compute, 'MxComputeDefault'))
	relevantOptions <- list(base::options()$mxOption$"Calculate Hessian", base::options()$mxOption$"Standard Errors",
													base::options()$mxOption$"Default optimizer")
	if("Calculate Hessian" %in%  names(model@options)){relevantOptions[[1]] <- model@options$"Calculate Hessian"}
	if("Standard Errors" %in%  names(model@options)){relevantOptions[[2]] <- model@options$"Standard Errors"}
	#If the options call for SEs and/or Hessian, there is no custom compute plan, and the Hessian will not be checked
	#every fit attempt, then computing SEs and/or Hessian can be put off until the MLE is obtained:
	SElater <- ifelse( (!checkHess && relevantOptions[[2]]=="Yes" && defaultComputePlan), TRUE, FALSE )
	Hesslater <- ifelse( (!checkHess && relevantOptions[[1]]=="Yes" && defaultComputePlan), TRUE, FALSE )
	if(SElater && !Hesslater){
		warning('the "Standard Errors" option is enabled and the "Calculate Hessian" option is disabled, which may result in poor-accuracy standard errors')
	}
	doIntervals <- ifelse ( (length(model@intervals) && intervals), TRUE, FALSE )
	lastNoError<-FALSE
	generalTolerance <- 1e-5 #used for hessian check and lowest min check
	gradientStepSize <- initialGradientStepSize
	tolerance <- initialTolerance
	gradientIterations <- initialGradientIterations  
	lastBestFitCount<-0 #number of consecutive improvements in fit
	stopflag <- FALSE #should the iterative optimization process stop?
	goodflag <- FALSE #is the best fit so far acceptable?
	numdone <- 0
	lowestminsofar<-Inf 
	finalfit<- NULL
	inits <- omxGetParameters(model)
	params <- inits
	if(is.na(maxMajorIter)){maxMajorIter <- max(1000, (3*length(inits)) + (10*length(model@constraints)))}
	parlbounds <- omxGetParameters(model=model,fetch="lbound")
	parlbounds[is.na(parlbounds)] <- -Inf
	parubounds <- omxGetParameters(model=model,fetch="ubound")
	parubounds[is.na(parubounds)] <- Inf
	
	
	#Begin main 'while' loop.
	while (!stopflag) {
		
		message(paste0('\nBegin fit attempt ', numdone+1, ' of at maximum ', extraTries +1, ' tries'))
		if(lastNoError && ("prev" %in% wtgcsv)){params <- omxGetParameters(fit)}
		
		
		if(lastBestFitCount == 0 && numdone > 0){ #if the last fit was not the best
			if(exists('bestfit') && ("best" %in% wtgcsv)){params <- bestfit.params} #if bestfit exists use this instead
			#sometimes, use initial start values instead:
			if(numdone %% 4 == 0 && ("initial" %in% wtgcsv)){params <- inits}
			model <- omxSetParameters(
				model, labels = names(params), 
				values=imxJiggle(params=params,lbounds=parlbounds,ubounds=parubounds,dsn=jitterDistrib,loc=loc,scale=scale)
			)
			if(finetuneGradient){
				gradientStepSize <- initialGradientStepSize
				tolerance <- initialTolerance
				gradientIterations<-initialGradientIterations
			}
		}#end if last fit not best section
		
		
		if(lastBestFitCount > 0){ #if the last fit was the best so far
			if(exists('bestfit')){
				if("best" %in% wtgcsv){params <- bestfit.params}
				model <- bestfit #<--Necessary?
			}
			if(defaultComputePlan==TRUE && finetuneGradient){
				if(lastBestFitCount == 2) gradientStepSize <- gradientStepSize *.1
				if(lastBestFitCount == 3) gradientStepSize <- gradientStepSize *10
				if(lastBestFitCount == 5) gradientStepSize <- gradientStepSize *.1
				if(lastBestFitCount  > 0) tolerance<-tolerance * .001 
				if(lastBestFitCount  > 0) gradientIterations<-gradientIterations+2
				if(lastBestFitCount > 2) model <- omxSetParameters(
					model, labels = names(bestfit.params), 
					values=imxJiggle(params=bestfit.params,lbounds=parlbounds,ubounds=parubounds,dsn=jitterDistrib,loc=loc,
													 scale=scale/10)
				)
			}
			else{
				model <- omxSetParameters(
					model, labels = names(bestfit.params), 
					values=imxJiggle(params=bestfit.params,lbounds=parlbounds,ubounds=parubounds,dsn=jitterDistrib,loc=loc,
													 scale=scale/ifelse(finetuneGradient,10,1))
				)
			}
		}#end if last fit was best section
		
		
		if(defaultComputePlan==TRUE){
			steps <- list(GD=mxComputeGradientDescent(
				verbose=verbose, gradientStepSize = gradientStepSize, 
				nudgeZeroStarts=FALSE,   gradientIterations = gradientIterations, tolerance=tolerance, 
				maxMajorIter=maxMajorIter))
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
		
		
		#If fit resulted in error:
		if( class(fit) == "try-error" || !is.finite(fit$output$minimum) || fit$output$status$status== -1){
			#^^^is.finite() returns FALSE for Inf, -Inf, NA, and NaN
			lastBestFitCount <- 0
			lastNoError<-FALSE
			message('\n Fit attempt generated errors') 
		}
		
		
		#If fit did NOT result in error:
		if(class(fit) != "try-error" && is.finite(fit$output$minimum) && fit$output$status$status != -1){
			lastNoError <- TRUE
			if(fit$output$minimum >= lowestminsofar){
				lastBestFitCount <- 0
				if(fit$output$minimum >= lowestminsofar + generalTolerance){
					message(paste0('\n Fit attempt worse than current best:  ',fit$output$minimum ,' vs ', lowestminsofar )) 
			}}
			#Current fit will become bestfit if (1) its fitvalue is strictly less than lowestminsofar, or
			#(2) its fitvalue is no greater than lowestminsofar (within tolerance) AND it satisfies the criteria for 
			#an acceptable result (i.e., goodflag gets set to TRUE):
			if(fit$output$minimum < lowestminsofar){ #<--If this is the best fit so far
				message(paste0('\n Lowest minimum so far:  ',fit$output$minimum) )
				lastBestFitCount<-lastBestFitCount+1 
				lowestminsofar <- fit$output$minimum
				bestfit <- fit
				bestfit.params <- omxGetParameters(bestfit)
			}
			if(fit$output$minimum <= lowestminsofar + generalTolerance){
				###########goodflag checks
				goodflag <- TRUE
				if( !(fit$output$status[[1]] %in% OKstatuscodes) ){
					goodflag <- FALSE
					message(paste0('\n OpenMx status code ', fit$output$status[[1]], ' not in list of acceptable status codes, ', OKstatuscodes))
				}
				if(fit$output$minimum > fit2beat) {
					message(paste0('\n Fit value of ', fit$output$minimum, ' greater than fit2beat of ', fit2beat))
					goodflag <- FALSE
				}
				if(checkHess==TRUE) {
					hessEigenval <- try(eigen(fit$output$calculatedHessian, symmetric = T, only.values = T)$values)
					if(class(hessEigenval)=='try-error') {
						message(paste0('\n Eigenvalues of Hessian could not be calculated'))
						goodflag <- FALSE
					}
					if(class(hessEigenval)!='try-error' && any(hessEigenval < 0)) {
						message(paste0('\n Not all eigenvalues of Hessian are greater than ', 0,': ', paste(hessEigenval,collapse=', ')))
						goodflag <- FALSE
				}}
				if(goodflag){ 
					bestfit <- fit
					bestfit.params <- omxGetParameters(bestfit)
				}
				stopflag <- goodflag && !exhaustive
			} #end goodflag checks
			
			if(iterationSummary){
				message(paste0("\n Attempt ",numdone," fit:  "))
				message(paste(names(params),": ", fit$output$estimate,"\n"))
				message(paste0("-2LL = ", fit$output$Minus2LogLikelihood))
			}
		} #end 'if fit did not result in error' section
		
		if(numdone > extraTries){
			message('\nRetry limit reached')
			stopflag <- TRUE
		}
	} #end while loop
	
	if(goodflag){
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
							maxMajorIter=maxMajorIter),
						constraintType=ifelse(relevantOptions[[3]] == 'NPSOL','none','ineq')))
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
				if (length(summary(finalfit)$npsolMessage) > 0){
					message('Warning messages generated from final fit for final fit for Hessian/SEs/CIs\n')
				}
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
	} #end 'if goodflag' section
	
	
	if(!goodflag){
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
								maxMajorIter=maxMajorIter),
							constraintType=ifelse(relevantOptions[[3]] == 'NPSOL','none','ineq')))
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
				if(bestfit$output$status$code==6 && !(6 %in% OKstatuscodes)){
					message('\nUncertain solution found - consider parameter validity, try again, increase extraTries, change inits, change model, or check data!\n')
				}
				if(iterationSummary==TRUE){
					message(paste(names(bestfit.params),": ", bestfit$output$estimate,"\n"))
					message(paste0("-2LL = ", bestfit$output$Minus2LogLikelihood))
				}
				bestfit <- THFrankenmodel(finalfit,bestfit,defaultComputePlan,Hesslater,SElater,doIntervals,checkHess)
			}}}
	
	
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
##' @param params Numeric vector of current free parameter values.
##' @param lbounds Numeric vector of lower bounds on parameters.
##' @param ubounds Numeric vector of upper bounds on parameters.
##' @param dsn Character string naming the family of distribution to be used.
##' @param loc Numeric vector of location parameters (medians).
##' @param scale Numeric vector of scale parameters (see source code).
##' @aliases
##' imxJiggle
imxJiggle <- function(params, lbounds, ubounds, dsn, loc, scale){
	if( !(dsn %in% c("rnorm","runif","rcauchy")) ){stop("unrecognized value for argument 'dsn'")}
	loc <- as.numeric(loc[1])
	scale <- as.numeric(scale[1])
	if(scale<0){stop("negative value for argument 'scale'")}
	n <- length(params)
	if(dsn=="rnorm"){
		out <- params * rnorm(n=n,mean=loc,sd=scale) + rnorm(n=n,mean=0,sd=scale)
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


#Wrapper function to imitate original implementation of mxTryHard()--attempts to find good start values:
mxTryHardOrig <- function(model, finetuneGradient=FALSE, maxMajorIter=NA, wtgcsv=c("prev","best"), ...){
	return(mxTryHard(model=model,finetuneGradient=finetuneGradient,
									 maxMajorIter=maxMajorIter,wtgcsv=wtgcsv,...))
}


#Wrapper function faithful to Charlie Driver's SSCT-oriented changes:
mxTryHardctsem <- function(model, initialGradientStepSize = .00001, initialGradientIterations = 1,
													initialTolerance=1e-12,	jitterDistrib="rnorm", ...){
	return(mxTryHard(model=model,initialGradientStepSize==initialGradientStepSize,
									 initialGradientIterations=initialGradientIterations,
									 initialTolerance=initialTolerance,jitterDistrib=jitterDistrib,...))
}


#Wrapper function that uses mxTryHard() to try to search a wide region of the parameter space:
mxTryHardWideSearch <- function(model, finetuneGradient=FALSE, jitterDistrib="rcauchy", exhaustive=TRUE,
	wtgcsv="prev", ...){
	return(mxTryHard(model=model,finetuneGradient==finetuneGradient,
									 jitterDistrib=jitterDistrib,
									 exhaustive=exhaustive,wtgcsv=wtgcsv,...))
}


#Wrapper function tailored toward ordinal-threshold analyses (not too sure about this function...): 
mxTryHardOrdinal <- function(model, greenOK = TRUE,	checkHess = FALSE, finetuneGradient=FALSE, exhaustive=TRUE,
	OKstatuscodes=c(0,1,5,6), wtgcsv=c("prev","best"), ...){
	return(mxTryHard(model=model,greenOK=greenOK,checkHess=checkHess,finetuneGradient=finetuneGradient,
									 exhaustive=exhaustive,OKstatuscodes=OKstatuscodes,wtgcsv=wtgcsv,...))
}
	