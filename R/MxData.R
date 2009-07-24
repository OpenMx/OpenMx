#
#   Copyright 2007-2009 The OpenMx Project
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
		matrix = "MxDataFrameOrMatrix",
		means  = "numeric",
		type   = "character",
		numObs = "numeric",
		name   = "character"))

setClassUnion("MxData", c("NULL", "MxNonNullData"))

setMethod("initialize", "MxNonNullData",
	function(.Object, matrix, means, type, numObs, name = "data") {
		.Object@matrix <- matrix
		.Object@means <- means
		.Object@type <- type
		.Object@numObs <- numObs
		.Object@name <- name
		return(.Object)
	}
)

omxDataTypes <- c("raw", "cov", "cor", "sscp")

mxData <- function(matrix, type, means = NA, numObs = NA) {
	if (length(means) == 1 && is.na(means)) means <- NA_real_
	if (missing(matrix) || !is(matrix, "MxDataFrameOrMatrix")) {
		stop("Matrix argument is neither a data frame nor a matrix")
	}
	if (missing(type) || (!is.character(type)) || (length(type) > 1) || 
		is.na(match(type, omxDataTypes))) {
		stop(paste("Type must be one of:", paste(omxDataTypes, collapse=" ")))
	}
	if (!is.vector(means) || !is.numeric(means)) {
		stop("Means argument must be of numeric vector type")
	}
	if (type != "raw" && is.na(numObs)) {
		stop("Number of observations must be specified for non-raw data")
	}
	if (type == "raw") {
		numObs <- nrow(matrix)
	} else if (type == 'cov' && (nrow(matrix) != ncol(matrix))) {
		stop("Matrix argument must be a square matrix for 'cov' data type")
	}
	numObs <- as.numeric(numObs)
	lapply(dimnames(matrix)[[2]], omxVerifyName)
	return(new("MxNonNullData", matrix, means, type, numObs))
}

displayMxData <- function(object) {
	cat("MxData", omxQuotes(object@name), '\n')
	cat("type :", omxQuotes(object@type), '\n')
	cat("numObs :", omxQuotes(object@numObs), '\n')
	cat("Matrix : \n") 
	print(object@matrix)
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
