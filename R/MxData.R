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
		indexVector = "integer",
		identicalDefVars = "integer",
		identicalMissingness = "integer",
		identicalRows = "integer",		
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
	lapply(dimnames(observed)[[2]], omxVerifyName, -1)
	means <- as.matrix(means)
	dim(means) <- c(1, length(means))
	return(new("MxNonNullData", observed, means, type, numObs))
}

convertDatasets <- function(model, defVars, modeloptions) {
	model@data <- sortRawData(model@data, defVars, model@name, modeloptions)
	model@data <- convertIntegerColumns(model@data)
	if (length(model@submodels) > 0) {
		model@submodels <- lapply(model@submodels, convertDatasets,
			defVars, modeloptions)
	}
	return(model)
}

sortRawData <- function(mxData, defVars, modelname, modeloptions) {
	if (is.null(mxData)) {
		return(mxData)
	}
	if (mxData@type != "raw") {
		return(mxData)	
	}
	observed <- mxData@observed
	nosort <- as.character(modeloptions[['No Sort Data']])
	fullname <- paste(modelname, 'data', sep = '.')
	components <- unlist(strsplit(fullname, omxSeparatorChar, fixed = TRUE))
	modelname <- components[[1]]	
	if ((length(observed) == 0) || (modelname %in% nosort)) {
		mxData@indexVector <- as.integer(NA)
		mxData@identicalDefVars <- as.integer(NA)
		mxData@identicalMissingness <- as.integer(NA)
		mxData@identicalRows <- as.integer(NA)			
	} else {
		observedNames <- colnames(observed)	
		if (length(defVars) > 0) {
			defKeys <- names(omxFilterDefinitionVariables(defVars, fullname))
			defKeys <- sapply(defKeys, function(x) {
				unlist(strsplit(x, omxSeparatorChar, fixed = TRUE))[[3]]
			})
			names(defKeys) <- NULL
			defKeys <- defKeys[defKeys %in% observedNames]
			defIndex <- match(defKeys, observedNames)
			otherIndex <- setdiff(1:length(observedNames), defIndex)
			sortkeys <- c(defIndex, otherIndex)
		} else {
			sortkeys <- c(1:length(observedNames))
			defKeys <- character()
		}
		sortvectors <- lapply(sortkeys, function(x) {observed[,x] })
		args <- c(sortvectors, 'na.last'=FALSE)
		indexVector <- do.call('order', args)
		sortdata <- observed[indexVector,,drop=FALSE]
		mxData@observed <- sortdata
		selectMissing <- is.na(sortdata)
		selectDefvars <- sortdata[, defKeys, drop=FALSE]
		threeVectors <- .Call("findIdenticalRowsData", sortdata, 
			selectMissing, selectDefvars, !any(selectMissing),
			length(selectDefvars) == 0, PACKAGE = "OpenMx")
		mxData@indexVector <- indexVector - 1L
		mxData@identicalRows <- threeVectors[[1]]
		mxData@identicalMissingness <- threeVectors[[2]]
		mxData@identicalDefVars <- threeVectors[[3]]
	}
	return(mxData)
}

calculateIdenticalDefVars <- function(sortdata, defKeys) {
	index <- 1
	retval <- integer()
	while(index <= nrow(sortdata)) {
		offset <- 1
		while((index + offset) <= nrow(sortdata) && 
			all(sortdata[index, defKeys] == sortdata[index + offset, defKeys])) {
			offset <- offset + 1
		}
		retval[[index]] <- as.integer(offset)
		if (offset > 1) { retval[(index + 1) : (index + offset - 1)] <- 0L }
		index <- index + offset
	}
	return(retval)
}

calculateIdenticalMissingness <- function(sortdata) {
	index <- 1
	retval <- integer()
	while(index <= nrow(sortdata)) {
		offset <- 1
		while((index + offset) <= nrow(sortdata) && 
			all(is.na(sortdata[index, ]) == is.na(sortdata[index + offset, ]))) {
			offset <- offset + 1
		}
		retval[[index]] <- as.integer(offset)
		if (offset > 1) { retval[(index + 1) : (index + offset - 1)] <- 0L }
		index <- index + offset
	}
	return(retval)
}

calculateIdenticalRows <- function(sortdata) {
	index <- 1
	retval <- integer()
	while(index <= nrow(sortdata)) {
		offset <- 1
		if (index + offset <= nrow(sortdata)) {
			leftside <- sortdata[index, ]
			rightside <- sortdata[index + offset, ]
			rownames(leftside) <- NULL
			rownames(rightside) <- NULL
			while(identical(leftside, rightside)) {
				offset <- offset + 1
				if (index + offset <= nrow(sortdata)) {
					rightside <- sortdata[index + offset, ]
					rownames(rightside) <- NULL
				} else {
					rightside <- NULL
				}
			}
		}
		retval[[index]] <- as.integer(offset)
		if (offset > 1) { retval[(index + 1) : (index + offset - 1)] <- 0L }
		index <- index + offset
	}
	return(retval)
}

convertIntegerColumns <- function(mxData) {
	if (is.null(mxData)) return(mxData)
	if (is.data.frame(mxData@observed)) {
		dimnames <- dimnames(mxData@observed)	
		mxData@observed <- data.frame(lapply(mxData@observed, convertIntegerSingleColumn))
		dimnames(mxData@observed) <- dimnames
	}
	mxData@numObs <- as.numeric(mxData@numObs)
	return(mxData)
}

convertIntegerSingleColumn <- function(column) {
	if(!is.factor(column) && is.integer(column)) {
		return(as.double(column))
	} else {
		return(column)
	}
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
