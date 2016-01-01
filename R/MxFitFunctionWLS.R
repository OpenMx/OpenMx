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

# TODO: Add conformance and dimname checking for weight matrix

#------------------------------------------------------------------------------
# Revision History
#  Fri 08 Mar 2013 15:35:29 Central Standard Time -- Michael Hunter copied file from 
#    OpenMx\branches\dependency-tracking\R\WLSObjective.R
#    I think Tim Brick wrote the original.
#    Michael Hunter edited the file to use fitfunctions
#  


#------------------------------------------------------------------------------

# **DONE**
setClass(Class = "MxFitFunctionWLS",
	representation = representation(weights = "MxCharOrNumber"),
	contains = "MxBaseFitFunction")

# **DONE**
setMethod("initialize", "MxFitFunctionWLS",
	function(.Object, weights, name = 'fitfunction') {
		.Object@name <- name
		.Object@weights <- weights
		.Object@vector <- FALSE
		return(.Object)
	}
)

# **DONE**
setMethod("qualifyNames", signature("MxFitFunctionWLS"), 
	function(.Object, modelname, namespace) {
		.Object@name <- imxIdentifier(modelname, .Object@name)
		.Object@weights <- imxConvertIdentifier(.Object@weights,
			modelname, namespace)
		return(.Object)
})

# **DONE**
setMethod("genericFitConvertEntities", "MxFitFunctionWLS", 
	function(.Object, flatModel, namespace, labelsData) {
		name <- .Object@name
		modelname <- imxReverseIdentifier(flatModel, .Object@name)[[1]]
		expectName <- paste(modelname, "expectation", sep=".")
		
		expectation <- flatModel@expectations[[expectName]]
		dataname <- expectation@data
		
		if (flatModel@datasets[[dataname]]@type != 'raw') {
			if (.Object@vector) {
				modelname <- getModelName(.Object)
				msg <- paste("The WLS fit function",
					"in model", omxQuotes(modelname),
					"does not have raw data")
				stop(msg, call.=FALSE)
			}
		}
		return(flatModel)
})

# **DONE**
setMethod("genericFitFunConvert", "MxFitFunctionWLS", 
	function(.Object, flatModel, model, labelsData, dependencies) {
		name <- .Object@name
		modelname <- imxReverseIdentifier(model, .Object@name)[[1]]
		expectName <- paste(modelname, "expectation", sep=".")
		if (expectName %in% names(flatModel@expectations)) {
			expectIndex <- imxLocateIndex(flatModel, expectName, name)
		} else {
			expectIndex <- as.integer(NA)
		}
		.Object@expectation <- expectIndex
		eobj <- flatModel@expectations[[1L+expectIndex]]
		if(!single.na(eobj@thresholds)){
			checkWLSIdentification(model, eobj)
		}
		return(.Object)
})

# **DONE**
setMethod("genericFitInitialMatrix", "MxFitFunctionWLS",
	function(.Object, flatModel) {
		flatFitFunction <- flatModel@fitfunctions[[.Object@name]]
		modelname <- imxReverseIdentifier(flatModel, flatFitFunction@name)[[1]]
		expectationName <- paste(modelname, "expectation", sep = ".")
		expectation <- flatModel@expectations[[expectationName]]
		if (is.null(expectation)) {
			msg <- paste("The WLS fit function",
			"has a missing expectation in the model",
			omxQuotes(modelname))
			stop(msg, call.=FALSE)
		}
		if (is.na(expectation@data)) {
			msg <- paste("The WLS fit function has",
			"an expectation function with no data in the model",
			omxQuotes(modelname))
			stop(msg, call.=FALSE)
		}
		mxDataObject <- flatModel@datasets[[expectation@data]]
		if (mxDataObject@type != 'raw') {
			msg <- paste("The dataset associated with the WLS expectation function", 
				"in model", omxQuotes(modelname), "is not raw data.")
			stop(msg, call.=FALSE)
		}
		cols <- ncol(mxDataObject@observed)
		return(matrix(as.double(NA), cols, cols))
})

setMethod("genericFitAddEntities", "MxFitFunctionWLS",
	function(.Object, job, flatJob, labelsData) {
		job <- mxOption(job, "Calculate Hessian", "No")
		job <- mxOption(job, "Standard Errors", "No")
		return(job)
})

checkWLSIdentification <- function(model, expect) {
		#TODO check the job's expectation
		# if it's expectation Normal
		#	if it's identified how we want, then do nothing
		#	else add constraints and change expectation
		#		might need to re-do expectationFunctionAddEntities and/or flattening
		# else if it's RAM
		#	if it's identified how we want, then do nothing
		#	else throw error
		# else if it's LISREL
		#	if it's identified how we want, then do nothing
		#	else throw error
		# else throw error that WLS doesn't work for that yet.

	type <- strsplit(is(expect)[1], "MxExpectation")[[1]][2]
	if(type=='RAM'){
		A <- expect@A
		S <- expect@S
		F <- expect@F
		M <- expect@M
		# ???
		# mxComputeOnce the RAM expectation?
		stop("Ordinal WLS with the RAM expectation is not yet implemented.  For now, use mxExpectationNormal.")
	}
	else if(type=='LISREL'){
		# ???
		# mxComputeOnce the LISREL expectation?
		stop("Ordinal WLS with the LISREL expectation is not yet implemented.  For now, use mxExpectationNormal.")
	}
	else if(type=='Normal'){
		#Note: below strategy does NOT catch cases where starting values
		# make the criteria true by accident.
		theCov <- mxEvalByName(expect@covariance, model, compute=TRUE) #Perhaps save time with cache?
		covDiagsOne <- all( ( (diag(theCov) - 1)^2 ) < 1e-7 )
		meanName <- expect@means
		if(!single.na(meanName)){
			theMeans <- mxEvalByName(expect@means, model, compute=TRUE)
			meanZeroNA <- all( ( (diag(theMeans) - 0)^2 ) < 1e-7 )
		} else {
			meanZeroNA <- TRUE
		}
		if(!covDiagsOne || !meanZeroNA){
			stop("Model not identified in the way required by WLS.")
		}
	}
	else{
		stop(paste("The model", omxQuotes(model@name), "has an expectation",
			"that is not RAM, LISREL, or Normal.  WLS does not know what to do."))
	}
}

# **DONE**
mxFitFunctionWLS <- function(weights="ULS") {
	if (is.na(weights)) weights <- as.character(NA)
	if (!(weights %in% c("ULS", "DWLS", "WLS"))){
		stop("weights must be either 'ULS', 'DLS' or 'WLS'")
		}
	return(new("MxFitFunctionWLS", weights))
}

# **DONE**
displayMxFitFunctionWLS <- function(fitfunction) {
	cat("MxFitFunctionWLS", omxQuotes(fitfunction@name), '\n')
	if (single.na(fitfunction@weights)) {
		cat("$weights : NA \n")
	} else {
		cat("$weights :", omxQuotes(fitfunction@weights), '\n')
	}
	if (length(fitfunction@result) == 0) {
		cat("$result: (not yet computed) ")
	} else {
		cat("$result:\n")
	}
	print(fitfunction@result)
	invisible(fitfunction)
}

# **DONE**
setMethod("print", "MxFitFunctionWLS", function(x, ...) { 
	displayMxFitFunctionWLS(x) 
})

setMethod("show", "MxFitFunctionWLS", function(object) { 
	displayMxFitFunctionWLS(object) 
})


imxWlsStandardErrors <- function(model){
	#TODO add safety check
	# Is it a WLS fit function
	# Does it have data of type=='acov'
	# Does the data have @fullWeight
	theParams <- omxGetParameters(model)
	d <- omxManifestModelByParameterJacobian(model)
	if(is.null(model$expectation) && (class(model$fitfunction) %in% "MxFitFunctionMultigroup") ){
		submNames <- sapply(strsplit(model$fitfunction$groups, ".", fixed=TRUE), "[", 1)
		sV <- list()
		sW <- list()
		for(amod in submNames){
			sV[[amod]] <- model[[amod]]$data$acov
			sW[[amod]] <- ginv(model[[amod]]$data$fullWeight)
		}
		V <- Matrix::bdiag(sV)
		W <- Matrix::bdiag(sW)
	} else {
		V <- model$data$acov #used weight matrix
		W <- MASS::ginv(model$data$fullWeight)
	}
	dvd <- solve( t(d) %*% V %*% d )
	nacov <- as.matrix(dvd %*% t(d) %*% V %*% W %*% V %*% d %*% dvd)
	wls.se <- matrix(sqrt(diag(nacov)), ncol=1)
	dimnames(nacov) <- list(names(theParams), names(theParams))
	rownames(wls.se) <- names(theParams)
	#SE is the standard errors
	#Cov is the analog of the Hessian for WLS
	return(list(SE=wls.se, Cov=nacov, Jac=d))
}



imxWlsChiSquare <- function(model, J=NA){
	samp.param <- mxGetExpected(model, 'vector')
	theParams <- omxGetParameters(model)
	if(is.null(model$expectation) && (class(model$fitfunction) %in% "MxFitFunctionMultigroup") ){
		submNames <- sapply(strsplit(model$fitfunction$groups, ".", fixed=TRUE), "[", 1)
		sW <- list()
		expd.param <- c()
		for(amod in submNames){
			cov <- model[[amod]]$data$observed
			mns <- model[[amod]]$data$means
			thr <- model[[amod]]$data$thresholds
			expd.param <- c(expd.param, cov[lower.tri(cov, TRUE)], mns[!is.na(mns)], thr[!is.na(thr)])
			sW[[amod]] <- MASS::ginv(model[[amod]]$data$fullWeight)
		}
		W <- Matrix::bdiag(sW)
	} else {
		cov <- model$data$observed
		mns <- model$data$means
		thr <- model$data$thresholds
		expd.param <- c(cov[lower.tri(cov, TRUE)], mns[!is.na(mns)], thr[!is.na(thr)])
		W <- MASS::ginv(model$data$fullWeight)
	}
	
	e <- samp.param - expd.param
	
	if(single.na(J)){
		jac <- omxManifestModelByParameterJacobian(model)
	} else {jac <- J}
	jacOC <- Null(jac)
	if(prod(dim(jacOC)) > 0){
		x2 <- t(e) %*% jacOC %*% ginv( as.matrix(t(jacOC) %*% W %*% jacOC) ) %*% t(jacOC) %*% e
	} else {x2 <- 0}
	df <- qr(jacOC)$rank
	return(list(Chi=x2, ChiDoF=df))
}



