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


setClass(Class = "SymmMatrix",
	representation = representation(),
	contains = "MxSymmetricMatrix")
	
setMethod("initialize", "SymmMatrix",
	function(.Object, name, values, free, labels, lbound, ubound, nrow, ncol, byrow) {
		if (nrow != ncol) {
			stop("Non-square matrix attempted for SymmMatrix constructor", call. = FALSE)
		}
		if (single.na(values)) {
			values <- 0
		}
		if (is.vector(values)) {
			len <- length(values)
			if (len == nrow * ncol || len == 1) {
				values <- matrix(values, nrow, ncol, byrow)
			} else if (len == nrow * (ncol + 1) / 2) {
				if(byrow) {
					tmp <- matrix(0, nrow, ncol)
					tmp[upper.tri(tmp, TRUE)] <- values
					tmp[lower.tri(tmp)] <- tmp[upper.tri(tmp)]
					values <- tmp
				} else {
					tmp <- matrix(0, nrow, ncol)
					tmp[lower.tri(tmp, TRUE)] <- values
					tmp[upper.tri(tmp)] <- tmp[lower.tri(tmp)]
					values <- tmp
				}				
			} else {
				stop(paste(
					"Illegal number of elements (", len,
					") for values matrix of SymmMatrix constructor", sep=""),
					call. = FALSE)
			}
		}
		if (is.vector(labels)) {
			len <- length(labels)
			if (len == nrow * ncol || len == 1) {
				labels <- matrix(labels, nrow, ncol, byrow)
			} else if (len == nrow * (ncol + 1) / 2) {
				if(byrow) {
					labels <- matrix(as.character(NA), nrow, ncol)
					tmp[upper.tri(tmp, TRUE)] <- labels
					tmp[lower.tri(tmp)] <- tmp[upper.tri(tmp)]
					labels <- tmp
				} else {
					tmp <- matrix(as.character(NA), nrow, ncol)
					tmp[lower.tri(tmp, TRUE)] <- labels
					tmp[upper.tri(tmp)] <- tmp[lower.tri(tmp)]
					labels <- tmp
				}				
			} else {
				stop(paste(
					"Illegal number of elements (", len,
					") for labels matrix of SymmMatrix constructor", sep=""),
					call. = FALSE)
			}
		}
		if (is.vector(free)) {
			len <- length(free)
			if (len == nrow * ncol || len == 1) {
				free <- matrix(free, nrow, ncol, byrow)
			} else if (len == nrow * (ncol + 1) / 2) {
				if(byrow) {
					free <- matrix(FALSE, nrow, ncol)
					tmp[upper.tri(tmp, TRUE)] <- free
					tmp[lower.tri(tmp)] <- tmp[upper.tri(tmp)]
					free <- tmp
				} else {
					free <- matrix(FALSE, nrow, ncol)
					tmp[lower.tri(tmp, TRUE)] <- free
					tmp[upper.tri(tmp)] <- tmp[lower.tri(tmp)]
					free <- tmp
				}				
			} else {
				stop(paste(
					"Illegal number of elements (", len,
					") for free matrix of SymmMatrix constructor", sep=""),
					call. = FALSE)
			}
		}
		if (is.vector(lbound)) {
			len <- length(lbound)
			if (len == nrow * ncol || len == 1) {
				lbound <- matrix(lbound, nrow, ncol, byrow)
			} else if (len == nrow * (ncol + 1) / 2) {
				if(byrow) {
					tmp <- matrix(as.numeric(NA), nrow, ncol)
					tmp[upper.tri(tmp, TRUE)] <- lbound
					tmp[lower.tri(tmp)] <- tmp[upper.tri(tmp)]
					lbound <- tmp
				} else {
					tmp <- matrix(as.numeric(NA), nrow, ncol)
					tmp[lower.tri(tmp, TRUE)] <- lbound
					tmp[upper.tri(tmp)] <- tmp[lower.tri(tmp)]
					lbound <- tmp
				}				
			} else {
				stop(paste(
					"Illegal number of elements (", len,
					") for lbound matrix of SymmMatrix constructor", sep=""),
					call. = FALSE)
			}
		}
		if (is.vector(ubound)) {
			len <- length(ubound)
			if (len == nrow * ncol || len == 1) {
				ubound <- matrix(ubound, nrow, ncol, byrow)
			} else if (len == nrow * (ncol + 1) / 2) {
				if(byrow) {
					tmp <- matrix(as.numeric(NA), nrow, ncol)
					tmp[upper.tri(tmp, TRUE)] <- ubound
					tmp[lower.tri(tmp)] <- tmp[upper.tri(tmp)]
					ubound <- tmp
				} else {
					tmp <- matrix(as.numeric(NA), nrow, ncol)
					tmp[lower.tri(tmp, TRUE)] <- ubound
					tmp[upper.tri(tmp)] <- tmp[lower.tri(tmp)]
					ubound <- tmp
				}				
			} else {
				stop(paste(
					"Illegal number of elements (", len,
					") for ubound matrix of SymmMatrix constructor", sep=""),
					call. = FALSE)
			}
		}
		retval <- callNextMethod(.Object, labels, values, free, lbound, ubound, name)
		return(retval)
	}
)

setMethod("omxVerifyMatrix", "SymmMatrix",
	function(.Object) {
		callNextMethod(.Object)
		verifySquare(.Object)
		values <- .Object@values
		free <- .Object@free
		labels <- .Object@labels
		lbound <- .Object@lbound
		ubound <- .Object@ubound
		if (!all(values == t(values))) {
			stop(paste("Values matrix of symmetric matrix", omxQuotes(.Object@name), 
				"is not symmetric!"), call. = FALSE)
		}
		if (!all(free == t(free))) {
			stop(paste("Free matrix of symmetric matrix", omxQuotes(.Object@name), 
				"is not symmetric!"), call. = FALSE)
		}
		if (!all(labels == t(labels), na.rm = TRUE) && all(is.na(labels) == is.na(t(labels)))) {
			stop(paste("Labels matrix of symmetric matrix", omxQuotes(.Object@name), 
				"is not symmetric!"), call. = FALSE)
		}
		if (!all(lbound == t(lbound), na.rm = TRUE) && all(is.na(lbound) == is.na(t(lbound)))) {
			stop(paste("Lbound matrix of symmetric matrix", omxQuotes(.Object@name), 
				"is not symmetric!"), call. = FALSE)
		}
		if (!all(ubound == t(ubound), na.rm = TRUE) && all(is.na(ubound) == is.na(t(ubound)))) {
			stop(paste("Ubound matrix of symmetric matrix", omxQuotes(.Object@name), 
				"is not symmetric!"), call. = FALSE)
		}
	}
)
