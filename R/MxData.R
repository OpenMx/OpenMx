#
#   Copyright 2007-2010 The OpenMx Project
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


setClassUnion("MxDataFrameOrMatrix", c("data.frame", "matrix"))

setClass(Class = "MxNonNullData",
	representation = representation(
		observed = "MxDataFrameOrMatrix",
		means  = "matrix",
		type   = "character",
		numObs = "numeric",
		name   = "character"))

setClassUnion("MxData", c("NULL", "MxNonNullData"))

setMethod("initialize", "MxNonNullData",
	function(.Object, observed, means, type, numObs, name = "data") {
		.Object@observed <- observed
		.Object@means <- means
		.Object@type <- type
		.Object@numObs <- numObs
		.Object@name <- name
		return(.Object)
	}
)

omxDataTypes <- c("raw", "cov", "cor", "sscp")

mxData <- function(observed, type, means = NA, numObs = NA) {
	if (length(means) == 1 && is.na(means)) means <- NA_real_
	if (missing(observed) || !is(observed, "MxDataFrameOrMatrix")) {
		stop("Observed argument is neither a data frame nor a matrix")
	}
	if (missing(type) || (!is.character(type)) || (length(type) > 1) || 
		is.na(match(type, omxDataTypes))) {
		stop(paste("Type must be one of:", paste(omxDataTypes, collapse=" ")))
	}
	if (type == "sscp") {
		stop(paste("'sscp' is not yet implemented."))
	}
	if (!is.vector(means) || !is.numeric(means)) {
		stop("Means argument must be of numeric vector type")
	}
	if (type != "raw" && is.na(numObs)) {
		stop("Number of observations must be specified for non-raw data")
	}
	if (type == "raw") {
		numObs <- nrow(observed)
	}
	if (type == "cov") {
		verifyCovarianceMatrix(observed)
	}
	numObs <- as.numeric(numObs)
	lapply(dimnames(observed)[[2]], omxVerifyName)
	means <- as.matrix(means)
	dim(means) <- c(1, length(means))
	return(new("MxNonNullData", observed, means, type, numObs))
}

convertIntegerSingleColumn <- function(column) {
	if(!is.factor(column) && is.integer(column)) {
		return(as.double(column))
	} else {
		return(column)
	}
}

convertIntegerColumns <- function(dframe) {
	output <- data.frame(lapply(dframe, convertIntegerSingleColumn))
	return(output)
}

checkNumericData <- function(data) {
	if(is.matrix(data@observed) && !is.double(data@observed)) {
		msg <- paste("The data object",
			omxQuotes(data@name), "contains an observed matrix that",
			"is not of type 'double'")
		stop(msg, call. = FALSE)
	}
}

verifyCovarianceMatrix <- function(covMatrix) {
	if(nrow(covMatrix) != ncol(covMatrix)) {
		msg <- paste("The observed covariance matrix",
			"is not a square matrix")
		stop(msg, call. = FALSE)
	}
	if (any(is.na(covMatrix))) {
		msg <- paste("The observed covariance matrix",
			"contains NA values")
		stop(msg, call. = FALSE)	
	}
	if (!all(covMatrix == t(covMatrix))) {
		msg <- paste("The observed covariance matrix",
			"is not a symmetric matrix")
		stop(msg, call. = FALSE)
	}
}

generateDataList <- function(model) {
	retval <- lapply(model@submodels, generateDataList)
	names(retval) <- NULL
	retval <- unlist(retval)
	if (is.null(retval)) retval <- list()
	retval[[model@name]] <- model@data
	return(retval)
}

undoDataShare <- function(model, dataList) {
	model@data <- dataList[[model@name]]
	model@submodels <- omxLapply(model@submodels, undoDataShare, dataList)
	return(model)
}


displayMxData <- function(object) {
	cat("MxData", omxQuotes(object@name), '\n')
	cat("type :", omxQuotes(object@type), '\n')
	cat("numObs :", omxQuotes(object@numObs), '\n')
	cat("Data Frame or Matrix : \n") 
	print(object@observed)
	if (length(object@means) == 1 && is.na(object@means)) {
		cat("Means : NA \n")
	} else {
		cat("Means : \n") 
		print(object@means)		
	}
	invisible(object)
}

setMethod("print", "MxNonNullData", function(x,...) { 
	displayMxData(x) 
})

setMethod("show", "MxNonNullData", function(object) { 
	displayMxData(object) 
})
