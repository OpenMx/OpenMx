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

#TODO:
#Need more input checking?  For instance, initialGradientIterations should be a positive integer, right?

runWithCounter <- function(model, count, silent, intervals=FALSE) {
	if(silent){
		plan <- model@compute
		plan <- mxComputeLoop(list(plan), i=count)
		fit <- mxRun(mxModel(model, plan), suppressWarnings = T,
								 silent=F, unsafe=T, intervals=intervals, beginMessage=FALSE)
		fit@compute <- fit@compute$steps[[1]]
		return(fit)
	}
	else{
		return(mxRun(model=model, suppressWarnings = T, unsafe=T, silent=T, intervals=intervals, beginMessage=T))
	}
}

mxTryHard <- function(
	model, extraTries = 10, greenOK = FALSE, loc = 1, 
	scale = 0.25,  initialGradientStepSize = imxAutoOptionValue("Gradient step size"), 
	initialGradientIterations = imxAutoOptionValue('Gradient iterations'),
	initialTolerance=as.numeric(mxOption(NULL,'Optimality tolerance')), 
	checkHess = TRUE, fit2beat = Inf, paste = TRUE,
	iterationSummary=FALSE, bestInitsOutput=TRUE, showInits=FALSE, verbose=0, intervals = FALSE,
	finetuneGradient=TRUE, jitterDistrib=c("runif","rnorm","rcauchy"), exhaustive=FALSE,
	maxMajorIter=3000, OKstatuscodes, wtgcsv=c("prev","best","initial"), silent=interactive()
){
	#Initialize stuff & check inputs:
	jitterDistrib <- match.arg(jitterDistrib)
	wtgcsv <- match.arg(wtgcsv,c("prev","best","initial"),several.ok=T)
	if(missing(OKstatuscodes)){OKstatuscodes <- as.integer(c(0,as.logical(greenOK[1])))}
	else if( !(0 %in% OKstatuscodes) ){OKstatuscodes <- c(OKstatuscodes,0)}
	#if( !("MxModel" %in% class(model)) ){stop("argument 'model' must be an object of class 'MxModel'")}
	if(initialTolerance<0){stop("value for argument 'initialTolerance' cannot be negative")}
  warnModelCreatedByOldVersion(model)
	if (omxHasDefaultComputePlan(model)) {
		model@compute <- NULL
	}
	lackOfConstraints <- !imxHasConstraint(model)
	hasThresholds <- imxHasThresholds(model)
	defaultComputePlan <- (is.null(model@compute) || is(model@compute, 'MxComputeDefault'))
	relevantOptions <- list(base::options()$mxOption$"Calculate Hessian", base::options()$mxOption$"Standard Errors",
													base::options()$mxOption$"Default optimizer", base::options()$mxOption$"Gradient algorithm")
	if("Calculate Hessian" %in%  names(model@options)){relevantOptions[[1]] <- model@options$"Calculate Hessian"}
	if("Standard Errors" %in%  names(model@options)){relevantOptions[[2]] <- model@options$"Standard Errors"}
	if("Gradient algorithm" %in%  names(model@options)){relevantOptions[[4]] <- model@options$"Gradient algorithm"}
	if(!lackOfConstraints){
		if(imxHasWLS(model)){relevantOptions[[2]] <- "No"}
		if(checkHess){
			message("Polite note from mxTryHard: Hessian not checked as model contains mxConstraints")
			checkHess <- FALSE
		}
	}
	ndgi <- ifelse(hasThresholds,3L,4L)
	ndgss <- ifelse(hasThresholds,1e-5,1e-7)
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
	previousLen <- 0L
	msg <- ""
	validcount <- 0
	errorcount <- 0
	fitvalAtStarts <- NA
	inits <- omxGetParameters(model)
	params <- inits
	if(is.na(maxMajorIter)){maxMajorIter <- max(1000, (3*length(inits)) + (10*length(model@constraints)))}
	parlbounds <- omxGetParameters(model=model,fetch="lbound")
	parlbounds[is.na(parlbounds)] <- -Inf
	parubounds <- omxGetParameters(model=model,fetch="ubound")
	parubounds[is.na(parubounds)] <- Inf
	
	
	#Get fit at start values:
	inputCompute <- model@compute
	model@compute <- mxComputeSequence(
		list(CO=mxComputeOnce(from="fitfunction", what="fit", .is.bestfit=TRUE),
		RE=mxComputeReportExpectation()))
	model@compute@.persist <- TRUE
	modelAtStartValues <- suppressWarnings(try(runWithCounter(model, 0, silent, F)))
	if(class(modelAtStartValues) != "try-error"){ 
		fitvalAtStarts <- modelAtStartValues@fitfunction@result[1]
		#If there are MxConstraints, we don't know if they're satisfied at the start values,
		#so we don't want to treat the fit at the start values as the lowest so far,
		#since an uphill step may be necessary to get to feasibility:
		if(is.finite(fitvalAtStarts) && lackOfConstraints){lowestminsofar <- fitvalAtStarts}
	}
	model@compute <- inputCompute
	rm(modelAtStartValues, inputCompute)
	

	#Begin main 'while' loop.
	while (!stopflag) {
		if(numdone==0){
			if(!silent){
				message("\nBeginning initial fit attempt")
			} else{
				msg <- "Beginning initial fit attempt"
				imxReportProgress(msg, previousLen)
				previousLen <- nchar(msg)
			}
		} else{
			if(!silent){
				message(paste0('\nBeginning fit attempt ', numdone, ' of at maximum ', extraTries, ' extra tries'))
			} else{
				msg <- paste0('Beginning fit attempt ', numdone, ' of at maximum ', extraTries, ' extra tries')
				imxReportProgress(msg, previousLen)
				previousLen <- nchar(msg)
			}
		}
		
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
				maxMajorIter=maxMajorIter, gradientAlgo=relevantOptions[[4]]))
			if(checkHess){steps <- c(steps,ND=mxComputeNumericDeriv(stepSize=ndgss,iterations=ndgi),SE=mxComputeStandardError())}
			model <- OpenMx::mxModel(
				model,
				mxComputeSequence(c( steps,RD=mxComputeReportDeriv(),RE=mxComputeReportExpectation() )))
		}
		#showInits=FALSE by default for mxTryHard() and its 4 specialized wrappers, and the extra printing that occurs when showInits=TRUE
		#is too much to summarize in one line.  Therefore, if the user has provided showInits=TRUE, then give him/her the extra printing 
		#requested notwithstanding the value of argument 'silent' (which by default is TRUE in an interactive session):
		if(showInits) {
			message('\nStarting values:  ')
			message(paste0(names(omxGetParameters(model)),' : ', omxGetParameters(model),'\n'))
		}
		fit <- suppressWarnings(try(runWithCounter(model, numdone, silent, intervals=F)))
		numdone <- numdone + 1
		
		
		#If fit resulted in error:
		if( class(fit) == "try-error" || !is.finite(fit@fitfunction@result[1]) || fit$output$status$status== -1){
			#^^^is.finite() returns FALSE for Inf, -Inf, NA, and NaN
			lastBestFitCount <- 0
			lastNoError<-FALSE
			errorcount <- errorcount + 1
			if(!silent){message('\n Fit attempt generated errors')}
		}
		
		
		#If fit did NOT result in error:
		if(class(fit) != "try-error" && is.finite(fit@fitfunction@result[1]) && fit$output$status$status != -1){
			lastNoError <- TRUE
			validcount <- validcount + 1
			if(fit@fitfunction@result[1] >= lowestminsofar){
				lastBestFitCount <- 0
				if(fit@fitfunction@result[1] >= lowestminsofar + generalTolerance){
					if(!silent){message(paste0('\n Fit attempt worse than current best:  ',fit@fitfunction@result[1] ,' vs ', lowestminsofar ))}
					else{
						msg <- paste0('Fit attempt ',numdone-1,', fit=',fit@fitfunction@result[1],', worse than previous best (',lowestminsofar,')')
						imxReportProgress(msg, previousLen)
						previousLen <- nchar(msg)
					}
				}}
			#Current fit will become bestfit if (1) its fitvalue is strictly less than lowestminsofar, or
			#(2) its fitvalue is no greater than lowestminsofar (within tolerance) AND it satisfies the criteria for 
			#an acceptable result (i.e., goodflag gets set to TRUE):
			if(fit@fitfunction@result[1] < lowestminsofar){ #<--If this is the best fit so far
				if(!silent){message(paste0('\n Lowest minimum so far:  ',fit@fitfunction@result[1]))}
				else{
					msg <- paste0('Fit attempt ',numdone-1,', fit=',fit@fitfunction@result[1],', new current best! (was ',lowestminsofar,')')
					imxReportProgress(msg, previousLen)
					previousLen <- nchar(msg)
				}
				lastBestFitCount<-lastBestFitCount+1 
				lowestminsofar <- fit@fitfunction@result[1]
				bestfit <- fit
				bestfit.params <- omxGetParameters(bestfit)
			}
			if(fit@fitfunction@result[1] <= lowestminsofar + generalTolerance){
				###########goodflag checks
				goodflag <- TRUE
				if( !(fit$output$status[[1]] %in% OKstatuscodes) ){
					goodflag <- FALSE
					if(!silent){
						message(paste(' OpenMx status code ', fit$output$status[[1]], ' not in list of acceptable status codes, ', 
													paste("(",paste(OKstatuscodes,collapse=","),")",sep=""), sep=""))
					}
				}
				if(fit@fitfunction@result[1] > fit2beat) {
					if(!silent){message(paste0(' Fit value of ', fit@fitfunction@result[1], ' greater than fit2beat of ', fit2beat))}
					goodflag <- FALSE
				}
				if(checkHess==TRUE) {
					fit@output["infoDefinite"] <- TRUE
					hessEigenval <- try(eigen(fit$output$calculatedHessian, symmetric = T, only.values = T)$values,silent=T)
					if(class(hessEigenval)=='try-error') {
						if(!silent){message(paste0(' Eigenvalues of Hessian could not be calculated'))}
						goodflag <- FALSE
						fit@output["infoDefinite"] <- FALSE
					}
					if(class(hessEigenval)!='try-error' && any(hessEigenval < 0)) {
						if(!silent){message(paste0(' Not all eigenvalues of the Hessian are positive: ', paste(hessEigenval,collapse=', ')))}
						goodflag <- FALSE
						fit@output["infoDefinite"] <- FALSE
					}}
				if(goodflag){ 
					bestfit <- fit
					bestfit.params <- omxGetParameters(bestfit)
				}
				stopflag <- goodflag && !exhaustive
			} #end goodflag checks
			
			#iterationSummary=FALSE by default for mxTryHard() and its 4 specialized wrappers, and the extra printing that occurs when 
			#iterationSummary=TRUE is too much to summarize in one line.  Therefore, if the user has provided iterationSummary=TRUE, then give him/her 
			#the extra printing requested notwithstanding the value of argument 'silent' (which by default is TRUE in an interactive session):
			if(iterationSummary){
				message(paste0("\n Attempt ",numdone-1," result:  "))
				message(paste(names(params),": ", fit$output$estimate,"\n"))
				message(paste0("fit value = ", fit@fitfunction@result[1]))
			}
		} #end 'if fit did not result in error' section
		
		if(numdone > extraTries){
			if(!silent){message('\nRetry limit reached')}
			stopflag <- TRUE
		}
	} #end while loop
	
	if(goodflag){
		if(!silent){message('\nSolution found\n')}
		if(any(Hesslater,SElater,doIntervals)){
			if(!silent){message("Final run, for Hessian and/or standard errors and/or confidence intervals\n")}
			else{
				msg <- 'Final run, for Hessian and/or standard errors and/or confidence intervals'
				imxReportProgress(msg, previousLen)
				previousLen <- nchar(msg)
			}
			finalfit <- bestfit
			if(defaultComputePlan){
				steps <- list()
				if(doIntervals){
					ciOpt <- mxComputeGradientDescent(
						nudgeZeroStarts=FALSE,gradientIterations=gradientIterations,
						tolerance=tolerance, maxMajorIter=maxMajorIter, gradientAlgo=relevantOptions[[4]])
					steps <- c(steps,CI=mxComputeConfidenceInterval(
						plan=ciOpt, constraintType=ciOpt$defaultCImethod))
				}
				if(Hesslater){
					steps <- c(steps,ND=mxComputeNumericDeriv(stepSize=ndgss,iterations=ndgi))
				} else {
					steps <- c(steps,ND=mxComputeNumericDeriv(knownHessian=bestfit$output$hessian,
																										checkGradient=FALSE,stepSize=ndgss,iterations=ndgi))
				}
				if(SElater){
					steps <- c(steps,SE=mxComputeStandardError(),HQ=mxComputeHessianQuality())
				}
				steps <- c(steps,RD=mxComputeReportDeriv())
				finalfit <- OpenMx::mxModel(finalfit,mxComputeSequence(steps=steps))
			}
			finalfit <- suppressWarnings(try(runWithCounter(finalfit, numdone, silent, intervals=doIntervals)))
			if(class(finalfit) == "try-error" || finalfit$output$status$status== -1) {
				if(!silent){message(' Errors during final fit for Hessian/SEs/CIs\n')}
			} else {
				if (length(summary(finalfit)$npsolMessage) > 0){
					if(!silent){message(' Warning messages generated from final run for Hessian/SEs/CIs\n')}
				}
			}
		}
		imxReportProgress("", previousLen)
		message(paste0("\n Solution found!  Final fit=", signif(bestfit@fitfunction@result[1],8), " (started at ", signif(fitvalAtStarts,8), ")  (" ,numdone, " attempt(s): ", validcount, " valid, ", errorcount," errors)\n"))
		if (length(summary(bestfit)$npsolMessage) > 0) {
			warning(summary(bestfit)$npsolMessage)
		}
		if(iterationSummary){
			message(paste(names(bestfit.params),": ", bestfit$output$estimate,"\n"))
			message(paste0("fit value = ", bestfit@fitfunction@result[1]))
		}
		bestfit <- THFrankenmodel(finalfit,bestfit,defaultComputePlan,Hesslater,SElater,doIntervals,checkHess,lackOfConstraints)
	} #end 'if goodflag' section
	
	
	if(!goodflag){
		if (exists("bestfit")) {
			if(any(Hesslater,SElater,doIntervals)){
				if(!silent){message("Computing Hessian and/or standard errors and/or confidence intervals from imperfect solution\n")}
				else{
					msg <- "Computing Hessian and/or standard errors and/or confidence intervals from imperfect solution"
					imxReportProgress(msg, previousLen)
					previousLen <- nchar(msg)
				}
				finalfit <- bestfit
				if(defaultComputePlan){
					steps <- list()
					if(doIntervals){
						ciOpt <- mxComputeGradientDescent(
							nudgeZeroStarts=FALSE,gradientIterations=gradientIterations,
							tolerance=tolerance, maxMajorIter=maxMajorIter, gradientAlgo=relevantOptions[[4]])
						steps <- c(steps,CI=mxComputeConfidenceInterval(
							plan=ciOpt, constraintType=ciOpt$defaultCImethod))
					}
					if(Hesslater){
						steps <- c(steps,ND=mxComputeNumericDeriv(stepSize=ndgss,iterations=ndgi))
					} else {
						steps <- c(steps,ND=mxComputeNumericDeriv(knownHessian=bestfit$output$hessian,
																											checkGradient=FALSE,stepSize=ndgss,iterations=ndgi))
					}
					if(SElater){
						steps <- c(steps,SE=mxComputeStandardError(),HQ=mxComputeHessianQuality())
					}
					steps <- c(steps,RD=mxComputeReportDeriv())
					finalfit <- OpenMx::mxModel(bestfit,mxComputeSequence(steps=steps))
				}
				finalfit <- suppressWarnings(try(runWithCounter(finalfit, numdone, silent, intervals=doIntervals)))
				if(class(finalfit) == "try-error" || finalfit$output$status$status== -1) {
					if(!silent){message('Errors occurred during final run for Hessian/SEs/CIs; returning best fit as-is\n')}
				}
			}
			imxReportProgress("", previousLen)
			message(paste0("\n Retry limit reached; Best fit=", signif(bestfit@fitfunction@result[1],8), " (started at ", signif(fitvalAtStarts,8), ")  (", numdone, " attempt(s): ", validcount, " valid, ", errorcount," errors)\n"))
			if (length(bestfit$output$status$statusMsg) > 0) { 
				warning(bestfit$output$status$statusMsg)
			}
			if(bestfit$output$status$code==6 && !(6 %in% OKstatuscodes)){
				if(!silent){message('\n Uncertain solution found - consider parameter validity, try again, increase extraTries, change inits, change model, or check data!\n')}
			}
			if(iterationSummary){
				message(paste(names(bestfit.params),": ", bestfit$output$estimate,"\n"))
				message(paste0("fit value = ", bestfit@fitfunction@result[1]))
			}
			bestfit <- THFrankenmodel(finalfit,bestfit,defaultComputePlan,Hesslater,SElater,doIntervals,checkHess,lackOfConstraints)
		}
	}
	
	
	if(bestInitsOutput && exists("bestfit")){
		bestfit.params <- omxGetParameters(bestfit)
		if(!silent){
			message(" Start values from best fit:")
			if(paste) message(paste(bestfit.params, sep=",", collapse = ",")) 
			if(!paste)  message(paste(names(bestfit.params),": ", bestfit.params,"\n"))
		}
	}
	
	if (!exists("bestfit")) {
		if(class(fit) == 'try-error') warning(fit[[length(fit)]])
		imxReportProgress("", previousLen)
		message('\n All fit attempts resulted in errors - check starting values or model specification\n')
		bestfit<-fit
	}
	
	if( defaultComputePlan && !("try-error" %in% class(bestfit)) ){bestfit@compute@.persist <- FALSE}
	
	return(bestfit)
}



imxJiggle <- function(params, lbounds, ubounds, dsn, loc, scale){
	if( !(dsn %in% c("rnorm","runif","rcauchy")) ){stop("unrecognized value for argument 'dsn'")}
	loc <- as.numeric(loc[1])
	scale <- as.numeric(scale[1])
	if(scale<0){stop("negative value for argument 'scale'")}
	lbounds[is.na(lbounds)] <- -Inf
	ubounds[is.na(ubounds)] <- Inf
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


mxJiggle <- function(model, classic=FALSE, dsn=c("runif","rnorm","rcauchy"), loc=1, scale=0.25){
  warnModelCreatedByOldVersion(model)
	dsn <- match.arg(dsn, c("runif","rnorm","rcauchy"))
	loc <- as.numeric(loc[1])
	scale <- as.numeric(scale[1])
	if(scale<0){stop("negative value for argument 'scale'")}
	params <- omxGetParameters(model=model,free=TRUE,fetch="values")
	labels <- names(params)
	lbounds <- omxGetParameters(model=model,free=TRUE,fetch="lbound")
	lbounds[is.na(lbounds)] <- -Inf
	ubounds <- omxGetParameters(model=model,free=TRUE,fetch="ubound")
	ubounds[is.na(ubounds)] <- Inf
	if(classic){
		out <- params + 0.1*(params+0.5)
		if(any(out<lbounds)){out[out<lbounds] <- lbounds[out<lbounds]}
		if(any(out>ubounds)){out[out>ubounds] <- ubounds[out>ubounds]}
		retval <- omxSetParameters(model=model,labels=labels,values=out)
	}
	else{
		retval <- omxSetParameters(
			model=model,labels=labels,
			values=imxJiggle(params=params,lbounds=lbounds,ubounds=ubounds,dsn=dsn,loc=loc,scale=scale)
			)
	}
	return(retval)
}



THFrankenmodel <- function(finalfit,bestfit,defaultComputePlan,Hesslater,SElater,doIntervals,checkHess,lackOfConstraints){
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
		bestfit@output$vcov <- finalfit@output$vcov
	}
	if(!lackOfConstraints){
		bestfit@output$constraintFunctionValues <- finalfit@output$constraintFunctionValues
		bestfit@output$constraintJacobian <- finalfit@output$constraintJacobian
		bestfit@output$constraintNames <- finalfit@output$constraintNames
		bestfit@output$constraintRows <- finalfit@output$constraintRows
		bestfit@output$constraintCols <- finalfit@output$constraintCols
	}
	bestfit@output$evaluations <- bestfit@output$evaluations + finalfit@output$evaluations
	bestfit@output$frontendTime <- bestfit@output$frontendTime + finalfit@output$frontendTime
	bestfit@output$backendTime <- bestfit@output$backendTime + finalfit@output$backendTime
	bestfit@output$independentTime <- bestfit@output$independentTime + finalfit@output$independentTime
	bestfit@output$wallTime <- bestfit@output$wallTime + finalfit@output$wallTime
	bestfit@output$cpuTime <- bestfit@output$cpuTime + finalfit@output$cpuTime
	needednames <- names(finalfit@output)[which(!(names(finalfit@output) %in% names(bestfit@output)))]
	#Whether or not output elements having to do with CIs, SEs, or Hessians should be returned is governed
	#by user-provided arguments to mxTryHard(); those elements are handled by code above:
	needednames <- needednames[
		!(needednames %in% c(
			"confidenceIntervals","confidenceIntervalCodes","calculatedHessian","hessian","standardErrors","infoDefinite","conditionNumber","vcov"))
	]
	bestfit@output[needednames] <- finalfit@output[needednames]
	bestfit@.modifiedSinceRun <- FALSE
	return(bestfit)
}


#Wrapper function to imitate original implementation of mxTryHard()--attempts to find good start values:
mxTryHardOrig <- function(model, finetuneGradient=FALSE, maxMajorIter=NA, wtgcsv=c("prev","best"), silent=FALSE, ...){
	return(mxTryHard(model=model,finetuneGradient=finetuneGradient,
									 maxMajorIter=maxMajorIter,wtgcsv=wtgcsv,silent=silent,...))
}


#Wrapper function faithful to Charlie Driver's SSCT-oriented changes:
mxTryHardctsem <- function(model, initialGradientStepSize = .00001, initialGradientIterations = 1,
													 initialTolerance=1e-12,	jitterDistrib="rnorm", ...){
	return(mxTryHard(model=model,initialGradientStepSize=initialGradientStepSize,
									 initialGradientIterations=initialGradientIterations,
									 initialTolerance=initialTolerance,jitterDistrib=jitterDistrib,...))
}


#Wrapper function that uses mxTryHard() to try to search a wide region of the parameter space:
mxTryHardWideSearch <- function(model, finetuneGradient=FALSE, jitterDistrib="rcauchy", exhaustive=TRUE,
																wtgcsv="prev", ...){
	return(mxTryHard(model=model,finetuneGradient=finetuneGradient,
									 jitterDistrib=jitterDistrib,
									 exhaustive=exhaustive,wtgcsv=wtgcsv,...))
}


#Wrapper function tailored toward ordinal-threshold analyses (not too sure about this function...): 
mxTryHardOrdinal <- function(model, greenOK = TRUE,	checkHess = FALSE, finetuneGradient=FALSE, exhaustive=TRUE,
														 OKstatuscodes=c(0,1,5,6), wtgcsv=c("prev","best"), ...){
	return(mxTryHard(model=model,greenOK=greenOK,checkHess=checkHess,finetuneGradient=finetuneGradient,
									 exhaustive=exhaustive,OKstatuscodes=OKstatuscodes,wtgcsv=wtgcsv,...))
}

