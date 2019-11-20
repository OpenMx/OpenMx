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


#------------------------------------------------------------------------------
# Author: Michael D. Hunter
# Date: 2016-01-30
# Filename: MxRobustSE.R
# Purpose: Write a function for robust standard errors
#------------------------------------------------------------------------------

##' imxRowGradients
##'
##' This is an internal function exported for those people who know
##' what they are doing.
##'
##' This function computes the gradient for each row of data.
##' The returned object is a matrix with the same number of rows as the data,
##' and the same number of columns as there are free parameters.
##'
##' @param model An OpenMx model object that has been run
##' @param robustSE Logical; are the row gradients being requested to calculate robust standard errors?
imxRowGradients <- function(model, robustSE=FALSE){
	if(is.null(model@output)){
		stop("The 'model' argument has no output.  Give me a model that has been run.")
	}
	#If 'model' contains submodels and uses the multigroup fitfunction, then we can be pretty sure what we're supposed to do to get
	#row gradients:
	if(length(model@submodels)){
		if(is(model@fitfunction,"MxFitFunctionMultigroup")){
			grads <- NULL
			paramLabels <- names(omxGetParameters(model))
			numParam <- length(paramLabels)
			contributingModelNames <- model@fitfunction$groups
			custom.compute <- mxComputeSequence(list(mxComputeNumericDeriv(checkGradient=FALSE, hessian=FALSE), mxComputeReportDeriv()))
			for(i in 1:length(model@submodels)){
				#Ignore submodels that don't contribute to the multigroup fitfunction:
				if(length(grep(pattern=model@submodels[[i]]$name,x=contributingModelNames))){
					currModel <- model@submodels[[i]]
					if(is.null(currModel@data)){
						if(robustSE){
							#This is a warning, not an error, because there are edge cases where the robust SEs could still be valid:
							warning(paste("submodel '",currModel@name,"' contributes to the multigroup fitfunction but contains no data; robust standard errors may be incorrect",sep=""))
						} 
						next
					}
					if(currModel@data$type!="raw"){
						if(robustSE){
							stop(paste("submodel '",currModel@name,"' contributes to the multigroup fitfunction but does not contain raw data, which is required for robust standard errors",sep=""))
						}
						next
					}
					if(robustSE && is(currModel@fitfunction, "MxFitFunctionWLS")){
						stop(paste("submodel '",currModel@name,"' contributes to the multigroup fitfunction but uses WLS fit; robust standard errors require ML fit",sep=""))	
					}
					if(length(currModel@submodels)){ #<--Possible TODO: handle this case with function recursion
						if(robustSE){
							warning(paste("submodel '",currModel@name,
														"' contains submodels of its own; support for submodels of submodels not implemented, so robust standard errors may be incorrect",sep=""))
						}
						else{
							warning(paste("submodel '",currModel@name,
														"' contains submodels of its own; support for submodels of submodels not implemented",sep=""))
						}
					}
					#By itself, a GREML model can't get robust SEs; you'd end up calculating the variance of the row derivatives for a sample of n=1 row.
					#But, if it's a submodel contributing to a multigroup fit (admittedly a corner case), then it's just another data row
					#(though I am not sure if the theory underlying the sandwich estimator still applies for restricted maximum likelihood):
					if(is(currModel@expectation, "MxExpectationGREML")){
						#There is assumed to be only one "row" with GREML expectation, even if the raw dataset isn't (yet) structured that way:
						currGrads <- matrix(0,nrow=1,ncol=numParam,dimnames=list(NULL,paramLabels))
						grun <- mxRun(mxModel(currModel, custom.compute))
						currGrads[1,names(grun$output$gradient)] <- grun$output$gradient
						grads <- rbind(grads,currGrads)
					}
					else{
						currData <- currModel@data$observed
						currGrads <- matrix(0,nrow=nrow(currData),ncol=numParam,dimnames=list(NULL,paramLabels))
						for(j in 1:nrow(currData)){
							grun <- mxRun(mxModel(currModel, custom.compute, mxData(currData[j,,drop=FALSE],"raw")), silent=as.logical((j-1)%%100))
							currGrads[j,names(grun$output$gradient)] <- grun$output$gradient
						}
						grads <- rbind(grads,currGrads)
					}
				}
			}
		}
		else{stop("to obtain gradients for data rows in submodels, please use an MxFitFunctionMultigroup in 'model'")}
	}
	else{ #i.e., if no submodels
		if(is.null(model@data)){
			stop("The 'model' argument must have data, or if multigroup, use an MxFitFunctionMultigroup")
		}
		if(model$data$type!='raw'){
			stop("The 'model' argument must have raw (not summary) data.")
		}
		nrows <- nrow(model$data$observed)
		data <- model@data@observed
		custom.compute <- mxComputeSequence(list(mxComputeNumericDeriv(checkGradient=FALSE, hessian=FALSE), mxComputeReportDeriv()))
		grads <- matrix(NA, nrows, length(coef(model)))
		gmodel <- model
		for(i in 1:nrows){
			gmodel <- mxModel(gmodel, custom.compute, mxData(data[i,,drop=FALSE], 'raw'))
			grun <- mxRun(gmodel, silent = as.logical((i-1)%%100), suppressWarnings = FALSE)
			grads[i,] <- grun$output$gradient #get gradient
		}
	}
	return(grads)
}


##' imxRobustSE
##'
##' This is an internal function exported for those people who know
##' what they are doing.
##'
##' This function computes robust standard errors via a sandwich estimator.
##' The "bread" of the sandwich is the numerically computed inverse Hessian
##' of the likelihood function.  This is what is typically used for standard
##' errors throughout OpenMx.  The "meat" of the sandwich is proportional to 
##' the covariance matrix of the numerically computed row derivatives of the 
##' likelihood function (i.e. row gradients).
##' 
##' When \code{details=FALSE}, only the standard errors are returned.
##' 
##' When \code{details=TRUE},
##' a list with five named elements is returned.  Element \code{SE} is the 
##' vector of standard errors that is also returned when \code{details=FALSE}.
##' Element \code{cov} is the full robust covariance matrix of the parameter 
##' estimates; the square root of the diagonal of \code{cov} gives the 
##' standard errors.  Element \code{bread} is the aforementioned 
##' "bread"--the naive (non-robust) covariance matrix of the parameter 
##' estimates.  Element \code{meat} is the aforementioned "meat," proportional
##' to the covariance matrix of the row gradients.  Element \code{TIC} 
##' is the model's Takeuchi Information Criterion, which is a generalization 
##' of AIC calculated from the "bread," the "meat," and the loglikelihood 
##' at the maximum-likelihood solution.
##' 
##' This function does not work correctly with multigroup models in which the 
##' groups themselves contain subgroups, or in which groups contain references
##' to objects in other groups.  This function also does not correctly handle
##' multilevel data.
##'
##' @param model An OpenMx model object that has been run.
##' @param details Logical. whether to return the full parameter
##' covariance matrix.
imxRobustSE <- function(model, details=FALSE){
	if(is(model@expectation, "MxExpectationGREML")){
		stop("robust standard errors cannot be calculated for a single-group model that uses GREML expectation")
	}
	if(!length(model@output$vcov)){
		stop("imxRobustSE() requires model to have a nonempty 'vcov' output slot (has the model been run?)")
	}
	if(imxHasWLS(model)){
		stop("'model' uses a WLS fitfunction; robust standard errors are automatically calculated for WLS models when they are run with the default compute plan")
	}
	if(!is(model@fitfunction, "MxFitFunctionML") && !is(model@fitfunction, "MxFitFunctionMultigroup")){
		warning(paste("imxRobustSE() requires a maximum-likelihood fit, but 'model' uses ",class(model@fitfunction),"; robust standard errors will only be correct if the fitfunction units are -2lnL",sep=""))
	}
	parnames <- dimnames(model@output$vcov)
	# if(!is.na(model@output$infoDefinite) && model@output$infoDefinite){
	# 	#solve() will fail if Hessian is computationally singular;
	# 	#chol2inv() will only fail if Hessian is exactly singular.
	# 	bread <- chol2inv(chol(model@output$hessian/2))
	# }
	# #An indefinite Hessian usually means some SEs will be NaN:
	# else{bread <- solve(model@output$hessian/2)}
	bread <- vcov(model)
	dimnames(bread) <- parnames
	#The row gradients are the slowest part, so only do them now that we know the bread is good:
	grads <- imxRowGradients(model, robustSE=TRUE)/-2
	meat <- nrow(grads)*var(grads)
	rm(grads) #<--Could be huge in Big Data contexts...
	dimnames(meat) <- parnames
	ret <- OpenMx::"%&%"(bread, meat)
	dimnames(ret) <- parnames
	TIC <- NA
	if(length(model@output$Minus2LogLikelihood)){TIC <- model@output$Minus2LogLikelihood + 2*sum(diag(meat%*%bread))}
	if(details){
		return(list(SE=sqrt(diag(ret)), cov=ret, bread=bread, meat=meat, TIC=TIC))
	} else {
		return(sqrt(diag(ret)))
	}
}

#robse <- imxRobustSE(thresholdModelrun)
#cbind(robse, prevSE)

