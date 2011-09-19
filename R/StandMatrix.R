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


setClass(Class = "StandMatrix",
	representation = representation(),
	contains = "MxMatrix")

setMethod("imxSymmetricMatrix", "StandMatrix",
	function(.Object) { return(TRUE) }
)

setMethod("imxSquareMatrix", "StandMatrix",
	function(.Object) { return(TRUE) }
)

populateStandTriangle <- function(input, n, default, byrow, strname) {
	len <- length(input)
	if (len == 1) {
		output <- matrix(default, n, n)
		output[row(output) != col(output)] <- input
	} else if (len == n * n) {
		output <- matrix(input, n, n, byrow)
	} else if (len == n * (n - 1) / 2) {
		if(byrow) {
			output <- matrix(default, n, n)
			output[upper.tri(output)] <- input
			output[lower.tri(output)] <- t(output)[lower.tri(output)]
		} else {
			output <- matrix(default, n, n)
			output[lower.tri(output)] <- input
			output[upper.tri(output)] <- t(output)[upper.tri(output)]
		}				
	} else {
		stop(paste(
			"Illegal number of elements (", len,
			") for ", strname, " matrix in standardized matrix constructor", sep="",
			deparse(width.cutoff = 400L, imxLocateFunction("mxMatrix"))),
			call. = FALSE)
	}
	return(output)
}
	
setMethod("imxCreateMatrix", "StandMatrix",
	function(.Object, labels, values, free, lbound, ubound, nrow, ncol, byrow, name, ...) {
		if (nrow != ncol) {
			stop("Non-square matrix attempted for standardized matrix constructor ",
			     deparse(width.cutoff = 400L, imxLocateFunction("mxMatrix")), 
                             call. = FALSE)
		}
		if (single.na(values)) {
			values <- 0
		}
		if (is.vector(values)) {
			values <- populateStandTriangle(values, nrow, 1, byrow, 'values') 
		}
		if (is.vector(labels)) {
			labels <- populateStandTriangle(labels, nrow, as.character(NA), byrow, 'labels')
		}
		if (is.vector(free)) {
			free <- populateStandTriangle(free, nrow, FALSE, byrow, 'free')
		}
		if (is.vector(lbound)) {
			lbound <- populateStandTriangle(lbound, nrow, as.numeric(NA), byrow, 'lbound')
		}
		if (is.vector(ubound)) {
			ubound <- populateStandTriangle(ubound, nrow, as.numeric(NA), byrow, 'ubound')
		}
		retval <- callNextMethod(.Object, labels, values, free, lbound, ubound, nrow, ncol, byrow, name, ...)
		return(retval)
	}
)

setMethod("imxVerifyMatrix", "StandMatrix",
	function(.Object) {
		callNextMethod(.Object)
		values <- .Object@values
		free <- .Object@free
		labels <- .Object@labels
		lbound <- .Object@lbound
		ubound <- .Object@ubound
		if (!all(values == t(values))) {
			stop(paste("Values matrix of standardized matrix", omxQuotes(.Object@name), 
				"is not symmetric!", 
				deparse(width.cutoff = 400L, imxLocateFunction("mxMatrix"))),
				call. = FALSE)
		}
		if (!all(diag(values) == 1)) {
			stop(paste("Values matrix of standardized matrix", omxQuotes(.Object@name), 
				"is not 1's along the diagonal!", 
				deparse(width.cutoff = 400L, imxLocateFunction("mxMatrix"))),
				call. = FALSE)
		}
		if (!all(free == t(free))) {
			stop(paste("Free matrix of standardized matrix", omxQuotes(.Object@name), 
				"is not symmetric!", 
				deparse(width.cutoff = 400L, imxLocateFunction("mxMatrix"))),
				call. = FALSE)
		}
		if (!all(diag(free) == FALSE)) {
			stop(paste("Free matrix of standardized matrix", omxQuotes(.Object@name), 
				"is not fixed along the diagonal!", 
				deparse(width.cutoff = 400L, imxLocateFunction("mxMatrix"))),
				call. = FALSE)
		}
		if (!all(labels == t(labels), na.rm = TRUE) && all(is.na(labels) == is.na(t(labels)))) {
			stop(paste("Labels matrix of standardized matrix", omxQuotes(.Object@name), 
				"is not symmetric!", 
				deparse(width.cutoff = 400L, imxLocateFunction("mxMatrix"))),
				call. = FALSE)
		}
		if (!all(is.na(diag(labels)))) {
			stop(paste("Labels matrix of standardized matrix", omxQuotes(.Object@name), 
				"is not NA along the diagonal!", 
				deparse(width.cutoff = 400L, imxLocateFunction("mxMatrix"))),
				call. = FALSE)
		}
		if (!all(lbound == t(lbound), na.rm = TRUE) && all(is.na(lbound) == is.na(t(lbound)))) {
			stop(paste("Lbound matrix of standardized matrix", omxQuotes(.Object@name), 
				"is not symmetric!", 
				deparse(width.cutoff = 400L, imxLocateFunction("mxMatrix"))),
				call. = FALSE)
		}
		if (!all(is.na(diag(lbound)))) {
			stop(paste("Lbound matrix of standardized matrix", omxQuotes(.Object@name), 
				"is not NA along the diagonal!", 
				deparse(width.cutoff = 400L, imxLocateFunction("mxMatrix"))),
				call. = FALSE)
		}
		if (!all(ubound == t(ubound), na.rm = TRUE) && all(is.na(ubound) == is.na(t(ubound)))) {
			stop(paste("Ubound matrix of standardized matrix", omxQuotes(.Object@name), 
				"is not symmetric!", 
				deparse(width.cutoff = 400L, imxLocateFunction("mxMatrix"))),
				call. = FALSE)
		}
		if (!all(is.na(diag(ubound)))) {
			stop(paste("Ubound matrix of standardized matrix", omxQuotes(.Object@name), 
				"is not NA along the diagonal!", 
				deparse(width.cutoff = 400L, imxLocateFunction("mxMatrix"))),
				call. = FALSE)
		}
	}
)
