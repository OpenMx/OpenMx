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


setClass(Class = "SdiagMatrix",
	representation = representation(),
	contains = "MxMatrix")

setMethod("omxSymmetricMatrix", "SdiagMatrix",
	function(.Object) { return(FALSE) }
)

setMethod("omxSquareMatrix", "SdiagMatrix",
	function(.Object) { return(TRUE) }
)
	
setMethod("initialize", "SdiagMatrix",
	function(.Object, name, values, free, labels, lbound, ubound, nrow, ncol, byrow) {
		if (nrow != ncol) {
			stop("Non-square matrix attempted for subdiagonal matrix constructor", call. = FALSE)
		}
		if (single.na(values)) {
			values <- 0
		}
		if (is.vector(values)) {
			len <- length(values)
			if (len == nrow * ncol) {
				values <- matrix(values, nrow, ncol, byrow)
			} else if (len == nrow * (ncol - 1) / 2 || len == 1) {
				if(byrow) {
					tmp <- matrix(0, nrow, ncol)
					tmp[upper.tri(tmp)] <- values
					tmp[lower.tri(tmp)] <- t(tmp)[lower.tri(tmp)]
					tmp[upper.tri(tmp, TRUE)] <- 0
					values <- tmp
				} else {
					tmp <- matrix(0, nrow, ncol)
					tmp[lower.tri(tmp)] <- values
					values <- tmp
				}			
			} else {
				stop(paste(
					"Illegal number of elements (", len,
					") for values matrix of subdiagonal matrix constructor", sep=""),
					call. = FALSE)
			}
		}
		if (is.vector(labels)) {
			len <- length(labels)
			if (len == nrow * ncol) {
				labels <- matrix(labels, nrow, ncol, byrow)
			} else if (len == nrow * (ncol - 1) / 2 || len == 1) {
				if(byrow) {
					tmp <- matrix(as.character(NA), nrow, ncol)
					tmp[upper.tri(tmp)] <- labels
					tmp[lower.tri(tmp)] <- t(tmp)[lower.tri(tmp)]
					tmp[upper.tri(tmp, TRUE)] <- as.character(NA)
					labels <- tmp
				} else {
					tmp <- matrix(as.character(NA), nrow, ncol)
					tmp[lower.tri(tmp)] <- labels
					labels <- tmp
				}				
			} else {
				stop(paste(
					"Illegal number of elements (", len,
					") for labels matrix of subdiagonal matrix constructor", sep=""),
					call. = FALSE)
			}
		}
		if (is.vector(free)) {
			len <- length(free)
			if (len == nrow * ncol) {
				free <- matrix(free, nrow, ncol, byrow)
			} else if (len == nrow * (ncol - 1) / 2 || len == 1) {
				if(byrow) {
					tmp <- matrix(FALSE, nrow, ncol)
					tmp[upper.tri(tmp)] <- free
					tmp[lower.tri(tmp)] <- t(tmp)[lower.tri(tmp)]
					tmp[upper.tri(tmp, TRUE)] <- FALSE
					free <- tmp
				} else {
					tmp <- matrix(FALSE, nrow, ncol)
					tmp[lower.tri(tmp)] <- free
					free <- tmp
				}
			} else {
				stop(paste(
					"Illegal number of elements (", len,
					") for free matrix of subdiagonal matrix constructor", sep=""),
					call. = FALSE)
			}
		}
		if (is.vector(lbound)) {
			len <- length(lbound)
			if (len == nrow * ncol) {
				lbound <- matrix(lbound, nrow, ncol, byrow)
			} else if (len == nrow * (ncol - 1) / 2 || len == 1) {
				if(byrow) {
					tmp <- matrix(as.numeric(NA), nrow, ncol)
					tmp[upper.tri(tmp)] <- lbound
					tmp[lower.tri(tmp)] <- t(tmp)[lower.tri(tmp)]
					tmp[upper.tri(tmp, TRUE)] <- as.numeric(NA)
					lbound <- tmp
				} else {
					tmp <- matrix(as.numeric(NA), nrow, ncol)
					tmp[lower.tri(tmp)] <- lbound
					lbound <- tmp
				}				
			} else {
				stop(paste(
					"Illegal number of elements (", len,
					") for lbound matrix of subdiagonal matrix constructor", sep=""),
					call. = FALSE)
			}
		}
		if (is.vector(ubound)) {
			len <- length(ubound)
			if (len == nrow * ncol) {
				ubound <- matrix(ubound, nrow, ncol, byrow)
			} else if (len == nrow * (ncol - 1) / 2 || len == 1) {
				if(byrow) {
					tmp <- matrix(as.numeric(NA), nrow, ncol)
					tmp[upper.tri(tmp)] <- ubound
					tmp[lower.tri(tmp)] <- t(tmp)[lower.tri(tmp)]
					tmp[upper.tri(tmp, TRUE)] <- as.numeric(NA)
					ubound <- tmp
				} else {
					tmp <- matrix(as.numeric(NA), nrow, ncol)
					tmp[lower.tri(tmp)] <- ubound
					ubound <- tmp
				}
			} else {
				stop(paste(
					"Illegal number of elements (", len,
					") for ubound matrix of subdiagonal matrix constructor", sep=""),
					call. = FALSE)
			}
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
