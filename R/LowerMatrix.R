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


setClass(Class = "LowerMatrix",
	representation = representation(),
	contains = "MxMatrix")

setMethod("omxSymmetricMatrix", "LowerMatrix",
	function(.Object) { return(FALSE) }
)

setMethod("omxSquareMatrix", "LowerMatrix",
	function(.Object) { return(TRUE) }
)
	
setMethod("initialize", "LowerMatrix",
	function(.Object, name, values, free, labels, lbound, ubound, nrow, ncol, byrow) {
		if (nrow != ncol) {
			stop("Non-square matrix attempted for lower matrix constructor", call. = FALSE)
		}
		if (single.na(values)) {
			values <- 0
		}
		if (is.vector(values)) {
			len <- length(values)
			if (len == nrow * ncol) {
				values <- matrix(values, nrow, ncol, byrow)
			} else if (len == nrow * (ncol + 1) / 2 || len == 1) {
				if(byrow) {
					tmp <- matrix(0, nrow, ncol)
					tmp[upper.tri(tmp, TRUE)] <- values
					tmp[lower.tri(tmp)] <- t(tmp)[lower.tri(tmp)]
					tmp[upper.tri(tmp)] <- 0
					values <- tmp
				} else {
					tmp <- matrix(0, nrow, ncol)
					tmp[lower.tri(tmp, TRUE)] <- values
					values <- tmp
				}			
			} else {
				stop(paste(
					"Illegal number of elements (", len,
					") for values matrix of lower matrix constructor", sep=""),
					call. = FALSE)
			}
		}
		if (is.vector(labels)) {
			len <- length(labels)
			if (len == nrow * ncol) {
				labels <- matrix(labels, nrow, ncol, byrow)
			} else if (len == nrow * (ncol + 1) / 2 || len == 1) {
				if(byrow) {
					tmp <- matrix(as.character(NA), nrow, ncol)
					tmp[upper.tri(tmp, TRUE)] <- labels
					tmp[lower.tri(tmp)] <- t(tmp)[lower.tri(tmp)]
					tmp[upper.tri(tmp)] <- as.character(NA)
					labels <- tmp
				} else {
					tmp <- matrix(as.character(NA), nrow, ncol)
					tmp[lower.tri(tmp, TRUE)] <- labels
					labels <- tmp
				}				
			} else {
				stop(paste(
					"Illegal number of elements (", len,
					") for labels matrix of lower matrix constructor", sep=""),
					call. = FALSE)
			}
		}
		if (is.vector(free)) {
			len <- length(free)
			if (len == nrow * ncol) {
				free <- matrix(free, nrow, ncol, byrow)
			} else if (len == nrow * (ncol + 1) / 2  || len == 1) {
				if(byrow) {
					tmp <- matrix(FALSE, nrow, ncol)
					tmp[upper.tri(tmp, TRUE)] <- free
					tmp[lower.tri(tmp)] <- t(tmp)[lower.tri(tmp)]
					tmp[upper.tri(tmp)] <- FALSE
					free <- tmp
				} else {
					tmp <- matrix(FALSE, nrow, ncol)
					tmp[lower.tri(tmp, TRUE)] <- free
					free <- tmp
				}
			} else {
				stop(paste(
					"Illegal number of elements (", len,
					") for free matrix of lower matrix constructor", sep=""),
					call. = FALSE)
			}
		}
		if (is.vector(lbound)) {
			len <- length(lbound)
			if (len == nrow * ncol) {
				lbound <- matrix(lbound, nrow, ncol, byrow)
			} else if (len == nrow * (ncol + 1) / 2 || len == 1) {
				if(byrow) {
					tmp <- matrix(as.numeric(NA), nrow, ncol)
					tmp[upper.tri(tmp, TRUE)] <- lbound
					tmp[lower.tri(tmp)] <- t(tmp)[lower.tri(tmp)]
					tmp[upper.tri(tmp)] <- as.numeric(NA)
					lbound <- tmp
				} else {
					tmp <- matrix(as.numeric(NA), nrow, ncol)
					tmp[lower.tri(tmp, TRUE)] <- lbound
					lbound <- tmp
				}				
			} else {
				stop(paste(
					"Illegal number of elements (", len,
					") for lbound matrix of lower matrix constructor", sep=""),
					call. = FALSE)
			}
		}
		if (is.vector(ubound)) {
			len <- length(ubound)
			if (len == nrow * ncol) {
				ubound <- matrix(ubound, nrow, ncol, byrow)
			} else if (len == nrow * (ncol + 1) / 2 || len == 1) {
				if(byrow) {
					tmp <- matrix(as.numeric(NA), nrow, ncol)
					tmp[upper.tri(tmp, TRUE)] <- ubound
					tmp[lower.tri(tmp)] <- t(tmp)[lower.tri(tmp)]
					tmp[upper.tri(tmp)] <- as.numeric(NA)
					ubound <- tmp
				} else {
					tmp <- matrix(as.numeric(NA), nrow, ncol)
					tmp[lower.tri(tmp, TRUE)] <- ubound
					ubound <- tmp
				}
			} else {
				stop(paste(
					"Illegal number of elements (", len,
					") for ubound matrix of lower matrix constructor", sep=""),
					call. = FALSE)
			}
		}
		retval <- callNextMethod(.Object, labels, values, free, lbound, ubound, name)
		return(retval)
	}
)

setMethod("omxVerifyMatrix", "LowerMatrix",
	function(.Object) {
		callNextMethod(.Object)
		values <- .Object@values
		free <- .Object@free
		labels <- .Object@labels
		lbound <- .Object@lbound
		ubound <- .Object@ubound
		if (!all(values[upper.tri(values)]  == 0)) {
			stop(paste("Upper triangle of values matrix in lower matrix", omxQuotes(.Object@name), 
				"is not all zeros!"), call. = FALSE)
		}
		if (!all(free[upper.tri(free)] == FALSE)) {
			stop(paste("Upper triangle of free matrix in lower matrix", omxQuotes(.Object@name), 
				"is not all fixed!"), call. = FALSE)
		}
		if (!all(is.na(labels[upper.tri(labels)]))) {
			stop(paste("Upper triangle of labels matrix in lower matrix", omxQuotes(.Object@name), 
				"is not all NAs!"), call. = FALSE)
		}
		if (!all(is.na(lbound[upper.tri(labels)]))) {
			stop(paste("Upper triangle of lbound matrix in lower matrix", omxQuotes(.Object@name), 
				"is not all NAs!"), call. = FALSE)
		}
		if (!all(is.na(ubound[upper.tri(labels)]))) {
			stop(paste("Upper triangle of ubound matrix in lower matrix", omxQuotes(.Object@name), 
				"is not all NAs!"), call. = FALSE)
		}
	}
)
