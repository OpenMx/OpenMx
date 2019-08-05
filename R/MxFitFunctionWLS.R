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
# Revision History
#  Fri 08 Mar 2013 15:35:29 Central Standard Time -- Michael Hunter copied file from 
#    OpenMx\branches\dependency-tracking\R\WLSObjective.R
#    I think Tim Brick wrote the original.
#    Michael Hunter edited the file to use fitfunctions
#  


#------------------------------------------------------------------------------

# **DONE**
setClass(Class = "MxFitFunctionWLS",
	contains = "MxBaseFitFunction",
	 representation = representation(
		 type = "character",
		 continuousType = "character",
		 fullWeight = "logical"
	))

# **DONE**
setMethod("initialize", "MxFitFunctionWLS",
	function(.Object, type, allContinuousMethod, fullWeight, name = 'fitfunction') {
		.Object@name <- name
		.Object@vector <- FALSE
		.Object@type <- type
		.Object@continuousType <- allContinuousMethod
		.Object@fullWeight <- fullWeight
		return(.Object)
	}
)

# **DONE**
setMethod("qualifyNames", signature("MxFitFunctionWLS"), 
	function(.Object, modelname, namespace) {
		.Object@name <- imxIdentifier(modelname, .Object@name)
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

# **DONE**
mxFitFunctionWLS <- function(type=c('WLS','DWLS','ULS'),
			     allContinuousMethod=c("cumulants", "marginals"),
			     fullWeight=TRUE) {
	type <- match.arg(type)
	allContinuousMethod <- match.arg(allContinuousMethod)
	return(new("MxFitFunctionWLS", type, allContinuousMethod, as.logical(fullWeight)))
}

# **DONE**
displayMxFitFunctionWLS <- function(fitfunction) {
	cat("MxFitFunctionWLS", omxQuotes(fitfunction@name), '\n')
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


# deprecated
# nocov start
imxWlsStandardErrors <- function(model){
	#TODO add safety check
	# Is it a WLS fit function
	# Does the data have @fullWeight
	isMultiGroupModel <- is.null(model$expectation) && (class(model$fitfunction) %in% "MxFitFunctionMultigroup")
	fwMsg <- "Terribly sorry, master, but you cannot compute standard errors without the full weight matrix."
	theParams <- omxGetParameters(model)
	if( isMultiGroupModel ){
		submNames <- sapply(strsplit(model$fitfunction$groups, ".", fixed=TRUE), "[", 1)
		sV <- list()
		sW <- list()
		sD <- c()
		for(amod in submNames){
			sV[[amod]] <- model[[amod]]$data$acov
			fullWeight <- model[[amod]]$data$fullWeight
			if(single.na(fullWeight)){stop(paste(fwMsg, '\nOffending model is', amod))}
			sW[[amod]] <- MASS::ginv(fullWeight)
			sD[[amod]] <- single.na(model[[amod]]$data$thresholds)
		}
		if( !(all(sD == TRUE) || all(sD == FALSE)) ){
			stop("I feel like I'm getting mixed signals.  You have some ordinal data, and some continuous data, and I'm not sure what to do.  Post this on the developer forums.")
		}
		d <- omxManifestModelByParameterJacobian(model, standardize=!any(sD))
		V <- Matrix::bdiag(sV)
		W <- Matrix::bdiag(sW)
	} else {
		d <- omxManifestModelByParameterJacobian(model, standardize=!single.na(model$data$thresholds))
		V <- model$data$acov #used weight matrix
		fullWeight <- model$data$fullWeight
		if(single.na(fullWeight)){stop(paste(fwMsg, '\nOffending model is', model@name))}
		W <- MASS::ginv(fullWeight)
	}
	dvd <- try(solve( t(d) %*% V %*% d ), silent=TRUE)
	if(class(dvd) != "try-error"){
		nacov <- as.matrix(dvd %*% t(d) %*% V %*% W %*% V %*% d %*% dvd)
		wls.se <- matrix(sqrt(diag(nacov)), ncol=1)
	} else {
		nacov <- matrix(NA, nrow=length(theParams), ncol=length(theParams))
		wls.se <- matrix(NA, nrow=length(theParams), ncol=1)
		warning("Standard error matrix is not invertible.\nIs your model not identified or in need of a constraint?\nCheck for identification with mxCheckIdentification.\nReturning model with all NA standard errors.")
	}
	dimnames(nacov) <- list(names(theParams), names(theParams))
	rownames(wls.se) <- names(theParams)
	#SE is the standard errors
	#Cov is the analog of the Hessian for WLS
	return(list(SE=wls.se, Cov=nacov, Jac=d))
}
# nocov end

# deprecated
# nocov start
imxWlsChiSquare <- function(model, J=NA){
	samp.param <- mxGetExpected(model, 'standVector')
	theParams <- omxGetParameters(model)
	numOrdinal <- 0
	numObs <- 0
	isMultiGroupModel <- is.null(model$expectation) && (class(model$fitfunction) %in% "MxFitFunctionMultigroup")
	fwMsg <- "Terribly sorry, master, but you cannot compute chi square without the full weight matrix."
	if( isMultiGroupModel ){
		submNames <- sapply(strsplit(model$fitfunction$groups, ".", fixed=TRUE), "[", 1)
		sV <- list()
		sW <- list()
		expd.param <- c()
		sD <- c()
		for(amod in submNames){
			cov <- model[[amod]]$data$observed
			mns <- model[[amod]]$data$means
			thr <- model[[amod]]$data$thresholds
			numObs <- numObs + model[[amod]]$data$numObs
			sD[[amod]] <- single.na(thr)
			if(!single.na(thr)){
				expd.param <- c(expd.param, .standardizeCovMeansThresholds(cov, mns, thr, !is.na(thr), vector=TRUE))
				numOrdinal <- numOrdinal + ncol(thr)
			} else {
				expd.param <- c(expd.param, cov[lower.tri(cov, TRUE)], mns[!is.na(mns)], thr[!is.na(thr)])
			}
			sV[[amod]] <- model[[amod]]$data$acov
			fullWeight <- model[[amod]]$data$fullWeight
			if(single.na(fullWeight)){stop(paste(fwMsg, '\nOffending model is', amod))}
			sW[[amod]] <- MASS::ginv(fullWeight)
		}
		V <- Matrix::bdiag(sV)
		W <- Matrix::bdiag(sW)
	} else {
		cov <- model$data$observed
		mns <- model$data$means
		thr <- model$data$thresholds
		numObs <- model$data$numObs
		sD <- single.na(thr)
		if(!single.na(thr)){
			expd.param <- .standardizeCovMeansThresholds(cov, mns, thr, !is.na(thr), vector=TRUE)
			numOrdinal <- numOrdinal + ncol(thr)
		} else {
			expd.param <- c(cov[lower.tri(cov, TRUE)], mns[!is.na(mns)], thr[!is.na(thr)])
		}
		V <- model$data$acov #used weight matrix
		fullWeight <- model$data$fullWeight
		if(single.na(fullWeight)){stop(paste(fwMsg, '\nOffending model is', model@name))}
		W <- MASS::ginv(fullWeight)
	}
	
	if( !(all(sD == TRUE) || all(sD == FALSE)) ){
		stop("I'm getting some mixed signals.  You have some ordinal data, and some continuous data, and I'm not sure what to do.  Post this on the developer forums.")
	}
	
	e <- samp.param - expd.param
	
	if(single.na(J)){
		jac <- omxManifestModelByParameterJacobian(model, standardize=!any(sD))
	} else {jac <- J}
	jacOC <- Null(jac)
	if(prod(dim(jacOC)) > 0){
		x2 <- t(e) %*% jacOC %*% ginv( as.matrix(t(jacOC) %*% W %*% jacOC) ) %*% t(jacOC) %*% e
	} else {x2 <- 0}
	df <- qr(jacOC)$rank
	
	dvd <- try(solve( t(jac) %*% V %*% jac ), silent=TRUE)
	if(class(dvd) != "try-error"){
		U <- V - V %*% jac %*% dvd %*% t(jac) %*% V
	} else {
		U <- matrix(NA, nrow=nrow(W), ncol=ncol(W))
	}
	UW <- as.matrix(U %*% W) # For multigroup case, convert Matrix type to matrix
	UW2 <- UW %*% UW # unclear if this should be UW^2 i.e. elementwise power
	trUW <- sum(diag(UW))
	madj <- trUW/df
	x2m <- as.numeric(model$fitfunction$result)/madj
	dstar <- round((trUW^2) / sum(diag(UW2)))
	mvadj <-  trUW^2/dstar
	x2mv <- as.numeric(model$fitfunction$result)/mvadj
	# N.B. x2mv is off by a factor of N where N is the total number of rows in all data sets for the ULS case.
	V <- as.matrix(V)
	I <- diag(1, nrow=nrow(V))
	x2mv <- x2mv*ifelse(all(V[V!=0] == I[V != 0]), 1/numObs, 1)
	return(list(Chi=x2, ChiDoF=df, ChiM=x2m, ChiMV=x2mv, mAdjust=madj, mvAdjust=mvadj, dstar=dstar))
}
# nocov end

approveWLSIntervals <- function(flatModel, modelName) {
	ff <- flatModel@fitfunctions[[ paste0(modelName, '.fitfunction') ]]
	if (is(ff, "MxFitFunctionWLS")) {
		if (ff@type != 'WLS') {
			stop(paste("Confidence intervals are not supported for DWLS or ULS. ",
				"Try mxSE or switch", omxQuotes(modelName), "to full WLS"))
		}
	} else if (is(ff, "MxFitFunctionMultigroup")) {
		for (g1 in flatModel@fitfunction$groups) {
			approveWLSIntervals(flatModel, g1)
		}
	}
}

#' Determine whether a dataset will have weights and summary statistics for the means if used with mxFitFunctionWLS
#'
#' Given either a data.frame or an mxData of type raw, this function determines whether \code{mxFitFunctionWLS}
#' will generate expectations for means.
#' 
#' All-continuous data processed using the "cumulants" method lack means, while
#' all continuous data processed with allContinuousMethod = "marginals" will have means.
#' 
#' When data are not all continuous, allContinuousMethod is ignored, and means are modelled.
#'
#' @param data the (currently raw) data being used in a \code{\link{mxFitFunctionWLS}} model.
#' @param allContinuousMethod the method used to process data when all columns are continuous.
#' @param verbose logical. Whether to report diagnostics.
#' @return - list describing the data.
#' @family Data Functions
#' @seealso - \code{\link{mxFitFunctionWLS}}, \code{\link{omxAugmentDataWithWLSSummary}}
#' @examples
#'
#' # ====================================
#' # = All continuous, data.frame input =
#' # ====================================
#'
#' tmp = mxDescribeDataWLS(mtcars, allContinuousMethod= "cumulants", verbose = TRUE)
#' tmp$hasMeans # FALSE - no means with cumulants
#' tmp = mxDescribeDataWLS(mtcars, allContinuousMethod= "marginals") 
#' tmp$hasMeans # TRUE we get means with marginals
#'
#' # ==========================
#' # = mxData object as input =
#' # ==========================
#' tmp = mxData(mtcars, type="raw")
#' mxDescribeDataWLS(tmp, allContinuousMethod= "cumulants", verbose = TRUE)$hasMeans # FALSE
#' mxDescribeDataWLS(tmp, allContinuousMethod= "marginals")$hasMeans  # TRUE
#'
#' # =======================================
#' # = One var is a factor: Means modelled =
#' # =======================================
#' tmp = mtcars
#' tmp$cyl = factor(tmp$cyl)
#' mxDescribeDataWLS(tmp, allContinuousMethod= "cumulants")$hasMeans # TRUE - always has means
#' mxDescribeDataWLS(tmp, allContinuousMethod= "marginals")$hasMeans # TRUE
#' 
mxDescribeDataWLS <- function(data, allContinuousMethod = c("cumulants", "marginals"), verbose=FALSE){
	allContinuousMethod = match.arg(allContinuousMethod)
	if(class(data) == "data.frame"){
		# all good
	} else if(class(data) == "MxDataStatic" && data$type == "raw"){
		data = data$observed
	}else{
		message("mxDescribeDataWLS currently only knows how to process dataframes and mxData of type = 'raw'.\n",
		"You offered up an object of class: ", omxQuotes(class(data)))
	}

	if(all(sapply(data, FUN= is.numeric))){
		if(verbose){ print("all continuous") }

		if(allContinuousMethod == "cumulants"){
			return(list(hasMeans = FALSE))
		} else {
			return(list(hasMeans = TRUE))
		}
	}else{
		# Data with any non-continuous vars have means under WLS
		return(list(hasMeans = TRUE))
	}
}


##' imxHasWLS
##'
##' This is an internal function exported for those people who know
##' what they are doing.  This function checks if a model uses a
##' fitfunction with WLS units.
##'
##' @param model model
imxHasWLS <- function(model){
	if(!is.null(model@output$fitUnits)){
		if(model@output$fitUnits=="r'Wr"){return(TRUE)}
		else{return(FALSE)}
	}
	if(is.null(model@fitfunction)){return(FALSE)}
	if(is(model@fitfunction, "MxFitFunctionWLS")){return(TRUE)}
	if(length(model@fitfunction$units) && model@fitfunction$units=="r'Wr"){return(TRUE)}
	if( is(model@fitfunction, "MXFitFunctionMultigroup") ){
		#Just in case the user provided 'modelName.fitfunction':
		submodnames <- unlist(lapply(strsplit(model@fitfunction@groups,"[.]"),function(x){x[1]}))
		for(i in 1:length(model@submodels)){
			if(model@submodels[[i]]@name %in% submodnames){
				probe <- imxHasWLS(model@submodels[[i]])
				if(probe){return(probe)}
			}
		}
	}
	return(FALSE)
}
