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

##' Build the model used for mxAutoStart
##'
##' @param model The MxModel for which starting values are desired
##' @param type The type of starting values to obtain, currently unweighted or diagonally weighted least squares, ULS or DWLS
##' @return
##' an MxModel that can be run to obtain starting values
##' @seealso \link{mxAutoStart}
omxBuildAutoStartModel <- function(model, type=c('ULS', 'DWLS')) {
	type <- match.barg(type)
	# Run the model through all the frontend processing to check for errors
	blah <- mxRun(model, silent=TRUE, suppressWarnings=FALSE, onlyFrontend=TRUE)
	blah <- NULL
	# If no errors found, continue with processing for autostart.
	if(is.null(model@fitfunction)){
		stop("I don't work with null fit functions.")
	}
	if(imxHasDefinitionVariable(model)){
		stop("Definition variables found. Automatic start values are not implemented for models with definition variables.\nNo plans have been made to add these soon.")
	}

	isMultiGroupModel <- is.null(model$expectation) && (class(model$fitfunction) %in% "MxFitFunctionMultigroup")

	if( isMultiGroupModel ){
		submNames <- sapply(strsplit(model$fitfunction$groups, ".", fixed=TRUE), "[", 1)
		wmodel <- model
		for(amod in submNames){
			wmodel[[amod]] <- mxModel(model[[amod]], autoStartDataHelper(model, subname=amod, type=type))
		}
		wmodel <- mxModel(wmodel, mxFitFunctionMultigroup(submNames))
	} else {
		wmodel <- mxModel(model, autoStartDataHelper(model, type=type))
	}
	wmodel <- mxOption(wmodel, "Calculate Hessian", "No")
	wmodel <- mxOption(wmodel, "Standard Errors", "No")
	wmodel
}

##' Automatically set starting values for an MxModel
##' 
##' @param model The MxModel for which starting values are desired
##' @param type The type of starting values to obtain, currently unweighted or diagonally weighted least squares, ULS or DWLS
##' 
##' @details
##' This function automatically picks very good starting values for many models (RAM, LISREL, Normal), including multiple group versions of these.
##' It works for models with algebras. Models of continuous, ordinal, and joint ordinal-continuous variables are also acceptable.
##' It works for models with covariance or raw data.
##' However, it does not currently work for models with definition variables, state space models, and item factor analysis models.
##' 
##' The method used to obtain new starting values is quite simple. The user's model is changed to an unweighted least squares (ULS) model. The ULS model is estimated and its final point estimates are returned as the new starting values. Optionally, diagonally weighted least squares (DWLS) can be used instead with the \code{type} argument.
##' 
##' Please note that ULS is sensitive to the scales of your variables. For example, if you have variables with means of 20 and variances of 0.001, then ULS will "weight" the means 20,000 times more than the variances and might result in zero variance estimates. Likewise if one variable has a variance of 20 and another has a variance of 0.001, the same problem may arise. To avoid this, make sure your variables are scaled accordingly. You could also use \code{type='DWLS'} to have the function use diagonally weighted least squares to obtain starting values.  Of course, using diagonally weighted least squares will take much much longer and will usually not provide better starting values than unweighted least squares.
##' 
##' Also note that if \code{model} contains a \link[=mxExpectationGREML]{GREML expectation}, argument \code{type} is ignored, and the function always uses a form of ULS.
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
mxAutoStart <- function(model, type=c('ULS', 'DWLS')){
  warnModelCreatedByOldVersion(model)
	if(is(model$expectation,"MxExpectationGREML")){
		return(autoStartGREML(model))
	}
	type <- match.barg(type)
	defaultComputePlan <- (is.null(model@compute) || is(model@compute, 'MxComputeDefault'))
	if(!defaultComputePlan){
		customPlan <- model@compute
		model@compute <- NULL
	}
	wmodel <- mxRun(omxBuildAutoStartModel(model, type), silent=TRUE)
	newparams <- coef(wmodel)
	oldparams <- coef(model)
	model <- omxSetParameters(model, values=newparams, labels=names(oldparams))
	if(!defaultComputePlan){
		model@compute <- customPlan
	}
	return(model)
}


#------------------------------------------------------------------------------

autoStartDataHelper <- function(model, subname=model@name, type){
	if(is.null(model[[subname]]@data)){
		stop(paste("Your model named", model[[subname]]@name, "doesn't have any data?  Sad."))
	}
	exps <- mxGetExpected(model, c('covariance', 'means', 'thresholds'), subname=subname)
	useVars <- dimnames(exps$covariance)[[1]]
	data <- model[[subname]]$data$observed
	origDataType <- model[[subname]]$data$type
	if(origDataType %in% 'cov'){
		data <- data[useVars,useVars]
		nrowData <- model[[subname]]$data$numObs
		meanData <- model[[subname]]$data$means
		mdata <- mxData(data, type=origDataType, numObs=nrowData, observedStats=list(cov=data))
		if (length(exps$means) > 0) {
			mdata$observedStats$means <- meanData
		}
	} else if (origDataType == 'raw') {
		data <- data[,useVars]
		# This conditional is for cases when the model has only 1 endogenous variable:
		if(!is.matrix(data) && !is.data.frame(data)){
			data <- as.matrix(data)
			colnames(data) <- useVars
		}
		if (type == 'ULS' && !any(sapply(data, is.ordered))) {
			# special case for ULS, all continuous
			os <- list(cov=cov(data, use='pair'))
			if (length(exps$means) > 0) os$means <- colMeans(data, na.rm=TRUE)
			return(list(mxData(data, type='raw', numObs=nrow(data), observedStats=os),
				mxFitFunctionWLS('ULS', fullWeight=FALSE)))
		}
		mdata <- mxData(data, 'raw')
	} else {
		stop(paste("mxAutoStart for",origDataType,"data is not implemented"))
	}
	list(mdata, mxFitFunctionWLS(type, ifelse(length(exps$means) > 0, 'marginals', 'cumulants'),
		type != 'ULS'))
}


autoStartGREML <- function(model){
	if(model$expectation$dataset.is.yX){
		y <- model$data$observed[,1]
		X <- model$data$observed[,-1]
		casesToDrop <- model$expectation$casesToDrop
	} else{
		dat <- mxGREMLDataHandler(
			data=model$data$observed,
			yvars=model$expectation$yvars,
			Xvars=model$expectation$Xvars,
			addOnes=model$expectation$addOnes,
			blockByPheno=model$expectation$blockByPheno,
			staggerZeroes=model$expectation$staggerZeroes
		)
		y <- dat$yX[,1]
		X <- dat$yX[,-1]
		casesToDrop <- dat$casesToDrop
		rm(dat)
	}
	olsresids <- lm(y~X+0)$residuals
	rm(X,y)
	Vdim <- nrow(mxEvalByName(model$expectation$V,model,T))
	if(length(olsresids) < Vdim){
		filt_mtx <- mxMatrix(type="Full",nrow=1,ncol=Vdim,free=F,values=1,name="filt")
		filt_mtx@values[casesToDrop] <- 0
		aff <- mxAlgebra( sum((S - omxSelectRowsAndCols(V,filt))%^%2), name="algfitfunc")
	} else{
		filt_mtx <- mxMatrix(type="Full",nrow=1,ncol=1,free=F,values=1,name="filt")
		aff <- mxAlgebra( sum((S-V)%^%2), name="algfitfunc")
	}
	tempmod <- mxModel(
		"tmp",
		model,
		aff,
		filt_mtx,
		mxMatrix(type="Symm",nrow=length(olsresids),free=F,values=outer(olsresids,olsresids),name="S",condenseSlots=T),
		mxAlgebraFromString(paste(model$name,model$expectation$V,sep="."),name="V"),
		mxFitFunctionAlgebra("algfitfunc")
	)
	rm(olsresids,aff,filt_mtx)
	tempmod <- mxOption(tempmod,"Calculate Hessian","No")
	tempmod <- mxOption(tempmod,"Standard Errors","No")
	tempmod <- mxRun(tempmod,silent=T)
	newparams <- coef(tempmod)
	rm(tempmod)
	oldparams <- coef(model)
	model <- omxSetParameters(model, values=newparams, labels=names(oldparams))
	return(model)
}
