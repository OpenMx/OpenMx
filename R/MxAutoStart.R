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


#------------------------------------------------------------------------------

##' Automatically set starting values for an MxModel
##' 
##' @param model The MxModel for which starting values are desired
##' 
##' @details
##' This function automatically picks very good starting values for many models (RAM, LISREL, Normal), including multiple group versions of these.
##' It works for models with algebras. Models of continuous, ordinal, and joint ordinal-continous variables are also acceptable.
##' It works for model with covariance or raw data.
##' However, it does not currently work for models with definition variables, state space models, and item factor analysis models.
##' 
##' The method used to obtain new starting values is quite simple. The user's model is changed to an unweighted least squares (ULS) model. The ULS model is estimated and its final point estimates are returned as the new starting values.
##' 
##' Please note that ULS is sensitive to the scales of your variables. For example, if you have variables with means of 20 and variances of 0.001, then ULS will "weight" the means 20,000 times more than the variances and might result in zero variance estimates. Likewise if one variable has a variance of 20 and another has a variance of 0.001, the same problem may arise. To avoid this, make sure your variables are scaled accordingly. You could also use diagonally weighted least squares to obtain your own starting values.
##' 
##' @return
##' an MxModel with new free parameter values
##' 
##' @examples
##' # Use the frontpage model with negative variances to show better
##' # starting values
##' library(OpenMx)
##' data(demoOneFactor)
##' 
##' latents  = c("G") # the latent factor
##' manifests = names(demoOneFactor) # manifest variables to be modeled
##' 
##' m1 <- mxModel("One Factor", type = "RAM", 
##' 	manifestVars = manifests, latentVars = latents, 
##' 	mxPath(from = latents, to = manifests),
##' 	mxPath(from = manifests, arrows = 2, values=-.2),
##' 	mxPath(from = latents, arrows = 2, free = FALSE, values = 1.0),
##' 	mxPath(from = "one", to = manifests),
##' 	mxData(demoOneFactor, type = "raw")
##' )
##' 
##' # Starting values imply negative variances!
##' mxGetExpected(m1, 'covariance')
##' 
##' # Use mxAutoStart to get much better starting values
##' m1s <- mxAutoStart(m1)
##' mxGetExpected(m1s, 'covariance')
mxAutoStart <- function(model){
	if(is.null(model@fitfunction)){
		stop("I don't work with null fit functions.")
	}
	if(imxHasDefinitionVariable(model)){
		stop("Definition variables found. Automatic start values are not implemented for models with definition variables.\nNo plans have been made to add these soon.")
	}
	
	isMultiGroupModel <- is.null(model$expectation) && (class(model$fitfunction) %in% "MxFitFunctionMultigroup")
	
	if( isMultiGroupModel ){
		submNames <- sapply(strsplit(model$fitfunction$groups, ".", fixed=TRUE), "[", 1)
		sD <- list()
		wmodel <- model
		for(amod in submNames){
			sD[[amod]] <- autoStartDataHelper(model, subname=amod)
			wmodel[[amod]] <- mxModel(model[[amod]], name=paste0('AutoStart', amod), sD[[amod]], mxFitFunctionWLS())
		}
		wmodel <- mxModel(wmodel, name='AutoStart', mxFitFunctionMultigroup(submNames))
	} else {
		mdata <- autoStartDataHelper(model)
		wmodel <- mxModel(model, name='AutoStart', mdata, mxFitFunctionWLS())
	}
	wmodel <- mxOption(wmodel, "Calculate Hessian", "No")
	wmodel <- mxOption(wmodel, "Standard Errors", "No")
	wmodel <- mxRun(wmodel, silent=TRUE)
	newparams <- coef(wmodel)
	oldparams <- coef(model)
	model <- omxSetParameters(model, values=newparams, labels=names(oldparams))
	return(model)
}


#------------------------------------------------------------------------------

autoStartDataHelper <- function(model, subname=model@name){
	if(is.null(model[[subname]]@data)){
		stop(paste("Your model named", model[[subname]]@name, "doesn't have any data?  Sad."))
	}
	exps <- mxGetExpected(model, c('covariance', 'means', 'thresholds'), subname=subname)
	wsize <- length(c(vech(exps$covariance), exps$means, exps$thresholds[!is.na(exps$thresholds)]))
	useVars <- dimnames(exps$covariance)[[1]]
	data <- model[[subname]]$data$observed[,useVars]
	hasOrdinal <- any(sapply(data, is.ordered))
	isCovData <- model[[subname]]$data$type %in% 'cov'
	I <- diag(1, wsize)
	if(isCovData){
		if (any(hasOrdinal)) {
			stop("Found ordinal data of type='cov'. I go crazy, crazy baby.")
		}
		covData <- data[useVars,]
		nrowData <- model[[subname]]$data$numObs
		meanData <- model[[subname]]$data$means
	} else {
		if(!hasOrdinal){
			covData <- cov(data, use='pair')
			nrowData <- nrow(data)
			meanData <- colMeans(data, na.rm=TRUE)
		} else {
			return(mxDataWLS(data, type="ULS", fullWeight=FALSE))
		}
	}
	mdata <- mxData(observed=I, type='acov', numObs=nrowData, 
			acov=I, fullWeight=I, means=meanData)
	mdata@observed <- covData
	return(mdata)
}
