#
#   Copyright 2007-20010 The OpenMx Project
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


setClass(Class = "SdiagMatrix",
	representation = representation(),
	contains = "MxMatrix")

setMethod("omxSymmetricMatrix", "SdiagMatrix",
	function(.Object) { return(FALSE) }
)

setMethod("omxSquareMatrix", "SdiagMatrix",
	function(.Object) { return(TRUE) }
)

populateSdiagTriangle <- function(input, n, default, byrow, strname) {
	len <- length(input)
	if (len == n * n) {
		output <- matrix(input, n, n, byrow)
	} else if (len == n * (n - 1) / 2 || len == 1) {
		if(byrow) {
			output <- matrix(default, n, n)
			output[upper.tri(output)] <- input
			output[lower.tri(output)] <- t(output)[lower.tri(output)]
			output[upper.tri(output, TRUE)] <- default
		} else {
			output <- matrix(default, n, n)
			output[lower.tri(output)] <- input
		}
	} else {
		stop(paste(
			"Illegal number of elements (", len,
			") for ", strname, " matrix of subdiagonal matrix constructor", sep=""),
			call. = FALSE)
	}
	return(output)
}
	
setMethod("initialize", "SdiagMatrix",
	function(.Object, name, values, free, labels, lbound, ubound, nrow, ncol, byrow) {
		if (nrow != ncol) {
			stop("Non-square matrix attempted for subdiagonal matrix constructor", call. = FALSE)
		}
		if (single.na(values)) {
			values <- 0
		}
		if (is.vector(values)) {
			values <- populateSdiagTriangle(values, nrow, 0, byrow, 'values')
		}
		if (is.vector(labels)) {
			labels <- populateSdiagTriangle(labels, nrow, as.character(NA), byrow, 'labels')
		}
		if (is.vector(free)) {
			free <- populateSdiagTriangle(free, nrow, FALSE, byrow, 'free')
		}
		if (is.vector(lbound)) {
			lbound <- populateSdiagTriangle(lbound, nrow, as.numeric(NA), byrow, 'lbound')
		}
		if (is.vector(ubound)) {
			ubound <- populateSdiagTriangle(ubound, nrow, as.numeric(NA), byrow, 'ubound')
		}
		retval <- callNextMethod(.Object, labels, values, free, lbound, ubound, name)
		return(retval)
	}
)

setMethod("omxVerifyMatrix", "SdiagMatrix",
	function(.Object) {
		callNextMethod(.Object)
		values <- .Object@values
		free <- .Object@free
		labels <- .Object@labels
		lbound <- .Object@lbound
		ubound <- .Object@ubound
		if (!all(values[upper.tri(values, TRUE)]  == 0)) {
			stop(paste("Upper triangle or diagonal of values matrix in subdiagonal matrix", omxQuotes(.Object@name), 
				"is not all zeros!"), call. = FALSE)
		}
		if (!all(free[upper.tri(free, TRUE)] == FALSE)) {
			stop(paste("Upper triangle or diagonal of free matrix in subdiagonal matrix", omxQuotes(.Object@name), 
				"is not all fixed!"), call. = FALSE)
		}
		if (!all(is.na(labels[upper.tri(labels, TRUE)]))) {
			stop(paste("Upper triangle or diagonal of labels matrix in subdiagonal matrix", omxQuotes(.Object@name), 
				"is not all NAs!"), call. = FALSE)
		}
		if (!all(is.na(lbound[upper.tri(labels, TRUE)]))) {
			stop(paste("Upper triangle or diagonal of lbound matrix in subdiagonal matrix", omxQuotes(.Object@name), 
				"is not all NAs!"), call. = FALSE)
		}
		if (!all(is.na(ubound[upper.tri(labels, TRUE)]))) {
			stop(paste("Upper triangle or diagonal of ubound matrix in subdiagonal matrix", omxQuotes(.Object@name), 
				"is not all NAs!"), call. = FALSE)
		}
	}
)
