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
		vector = "numeric",
		type   = "character",
		numObs = "numeric",
		name   = "character"))

setClassUnion("MxData", c("NULL", "MxNonNullData"))

setMethod("initialize", "MxNonNullData",
	function(.Object, matrix, vector, type, numObs, name = "data") {
		.Object@matrix <- matrix
		.Object@vector <- vector
		.Object@type <- type
		.Object@numObs <- numObs
		.Object@name <- name
		return(.Object)
	}
)

omxDataTypes <- c("raw", "cov", "cor", "sscp")

mxData <- function(matrix, type, vector = NA, numObs = NA) {
	if(is.na(vector)) vector <- NA_real_
	if(missing(matrix) || !is(matrix, "MxDataFrameOrMatrix")) {
		stop("Matrix argument is neither a data frame nor a matrix")
	}
	if(missing(type) || (!is.character(type)) || (length(type) > 1) || 
		is.na(match(type, omxDataTypes))) {
		stop(paste("Type must be one of:", paste(omxDataTypes, collapse=" ")))
	}
	if(!is.vector(vector) || !is.numeric(vector)) {
		stop("Vector argument must be of numeric vector type")
	}
	if(type != "raw" && is.na(numObs)) {
		stop("Number of observations must be specified for non-raw data")
	}
	if(type == "raw") {
		numObs <- nrow(matrix)
	}
	lapply(dimnames(matrix)[[2]], omxVerifyName)
	return(new("MxNonNullData", matrix, vector, type, numObs))
}

displayMxData <- function(object) {
	cat("MxData", omxQuotes(object@name), '\n')
	cat("type :", omxQuotes(object@type), '\n')
	cat("numObs :", omxQuotes(object@numObs), '\n')
	cat("Matrix : \n") 
	print(object@matrix)
	if (is.na(object@vector)) {
		cat("Vector : NA \n")
	} else {
		cat("Vector : \n") 
		print(object@vector)		
	}
	invisible(object)
}

setMethod("print", "MxNonNullData", function(x,...) { 
	displayMxData(x) 
})

setMethod("show", "MxNonNullData", function(object) { 
	displayMxData(object) 
})
