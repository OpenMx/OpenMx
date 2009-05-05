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
	function(.Object, matrix, vector, type, numObs, name) {
		.Object@matrix <- matrix
		.Object@vector <- vector
		.Object@type <- type
		.Object@numObs <- numObs
		.Object@name <- name
		return(.Object)
	}
)

omxDataTypes <- c("raw", "cov", "cor", "sscp")

mxData <- function(matrix, type, vector = NA, numObs = NA, name = NA) {
	if(is.na(name)) name <- omxUntitledName()
	omxVerifyName(name)
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
	if (typeof(name) != "character") {
		stop("Name argument is not a string (the name of the objective function)")
	}
	return(new("MxNonNullData", matrix, vector, type, numObs, name))
}
