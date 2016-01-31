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
imxRowGradients <- function(model){
	# TODO something with multigroup models
	if(is.null(model@data)){
		stop("The 'model' argument must have data.")
	}
	if(model$data$type!='raw'){
		stop("The 'model' argument must have raw (not summary) data.")
	}
	nrows <- nrow(model$data$observed)
	if(length(model$data$indexVector) == nrows){ #put data back in unsorted order
		model@data@observed <- model$data$observed[order(model$data$indexVector), ]
	}
	if(is.null(model@output)){
		stop("The 'model' argument has no output.  Give me a model that has been run.")
	}
	data <- model@data@observed
	custom.compute <- mxComputeSequence(list(mxComputeNumericDeriv(checkGradient=FALSE, hessian=FALSE), mxComputeReportDeriv()))
	grads <- matrix(NA, nrows, length(coef(model)))
	gmodel <- model
	for(i in 1:nrows){
		gmodel <- mxModel(gmodel, custom.compute, mxData(data[i,,drop=FALSE], 'raw'))
		grun <- mxRun(gmodel, silent = as.logical((i-1)%%100), suppressWarnings = FALSE)
		grads[i,] <- grun$output$gradient #get gradient
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
##' errors throughout OpenMx.  The "meat" of the sandwich is the covariance
##' matrix of the numerically computed row derivatives of the likelihood function
##' (i.e. row gradients).
##'
##' @param model An OpenMx model object that has been run
imxRobustSE <- function(model){
	grads <- imxRowGradients(model)/-2
	hess <- model@output$hessian/2
	ret <- OpenMx::"%&%"(solve(hess), nrow(grads)*var(grads))
	return(sqrt(diag(ret)))
}

#robse <- imxRobustSE(thresholdModelrun)
#cbind(robse, prevSE)

