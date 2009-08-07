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


setClass(Class = "StandMatrix",
	representation = representation(),
	contains = "MxMatrix")

setMethod("omxSymmetricMatrix", "StandMatrix",
	function(.Object) { return(TRUE) }
)

setMethod("omxSquareMatrix", "StandMatrix",
	function(.Object) { return(TRUE) }
)
	
setMethod("initialize", "StandMatrix",
	function(.Object, name, values, free, labels, lbound, ubound, nrow, ncol, byrow) {
		if (nrow != ncol) {
			stop("Non-square matrix attempted for standardized matrix constructor", call. = FALSE)
		}
		if (single.na(values)) {
			values <- 0
		}
		if (is.vector(values)) {
			len <- length(values)
			if (len == 1) {
				tmp <- diag(1, nrow, ncol)
				tmp[row(tmp) != col(tmp)] <- values
				values <- tmp
			} else if (len == nrow * ncol) {
				values <- matrix(values, nrow, ncol, byrow)
			} else if (len == nrow * (ncol - 1) / 2) {
				if(byrow) {
					tmp <- matrix(1, nrow, ncol)
					tmp[upper.tri(tmp)] <- values
					tmp[lower.tri(tmp)] <- t(tmp)[lower.tri(tmp)]
					values <- tmp
				} else {
					tmp <- matrix(1, nrow, ncol)
					tmp[lower.tri(tmp)] <- values
					tmp[upper.tri(tmp)] <- t(tmp)[upper.tri(tmp)]
					values <- tmp
				}				
			} else {
				stop(paste(
					"Illegal number of elements (", len,
					") for values matrix in standardized matrix constructor", sep=""),
					call. = FALSE)
			}
		}
		if (is.vector(labels)) {
			len <- length(labels)
			if (len == 1) {
				tmp <- diag(as.character(NA), nrow, ncol)
				tmp[row(tmp) != col(tmp)] <- labels
				labels <- tmp
			} else if (len == nrow * ncol) {
				labels <- matrix(labels, nrow, ncol, byrow)
			} else if (len == nrow * (ncol - 1) / 2) {
				if(byrow) {
					tmp <- matrix(as.character(NA), nrow, ncol)
					tmp[upper.tri(tmp)] <- labels
					tmp[lower.tri(tmp)] <- t(tmp)[lower.tri(tmp)]
					labels <- tmp
				} else {
					tmp <- matrix(as.character(NA), nrow, ncol)
					tmp[lower.tri(tmp)] <- labels
					tmp[upper.tri(tmp)] <- t(tmp)[upper.tri(tmp)]
					labels <- tmp
				}				
			} else {
				stop(paste(
					"Illegal number of elements (", len,
					") for labels matrix in standardized matrix constructor", sep=""),
					call. = FALSE)
			}
		}
		if (is.vector(free)) {
			len <- length(free)
			if (len == 1) {
				tmp <- diag(FALSE, nrow, ncol)
				tmp[row(tmp) != col(tmp)] <- free
				free <- tmp
			} else if (len == nrow * ncol) {
				free <- matrix(free, nrow, ncol, byrow)
			} else if (len == nrow * (ncol - 1) / 2) {
				if(byrow) {
					tmp <- matrix(FALSE, nrow, ncol)
					tmp[upper.tri(tmp)] <- free
					tmp[lower.tri(tmp)] <- t(tmp)[lower.tri(tmp)]
					free <- tmp
				} else {
					tmp <- matrix(FALSE, nrow, ncol)
					tmp[lower.tri(tmp)] <- free
					tmp[upper.tri(tmp)] <- t(tmp)[upper.tri(tmp)]
					free <- tmp
				}				
			} else {
				stop(paste(
					"Illegal number of elements (", len,
					") for free matrix in standardized matrix constructor", sep=""),
					call. = FALSE)
			}
		}
		if (is.vector(lbound)) {
			len <- length(lbound)
			if (len == 1) {
				tmp <- diag(as.numeric(NA), nrow, ncol)
				tmp[row(tmp) != col(tmp)] <- lbound
				lbound <- tmp
			} else if (len == nrow * ncol) {
				lbound <- matrix(lbound, nrow, ncol, byrow)
			} else if (len == nrow * (ncol - 1) / 2) {
				if(byrow) {
					tmp <- matrix(as.numeric(NA), nrow, ncol)
					tmp[upper.tri(tmp)] <- lbound
					tmp[lower.tri(tmp)] <- t(tmp)[lower.tri(tmp)]
					lbound <- tmp
				} else {
					tmp <- matrix(as.numeric(NA), nrow, ncol)
					tmp[lower.tri(tmp)] <- lbound
					tmp[upper.tri(tmp)] <- t(tmp)[upper.tri(tmp)]
					lbound <- tmp
				}				
			} else {
				stop(paste(
					"Illegal number of elements (", len,
					") for lbound matrix in standardized matrix constructor", sep=""),
					call. = FALSE)
			}
		}
		if (is.vector(ubound)) {
			len <- length(ubound)
			if (len == 1) {
				tmp <- diag(as.numeric(NA), nrow, ncol)
				tmp[row(tmp) != col(tmp)] <- ubound
				ubound <- tmp
			} else if (len == nrow * ncol) {
				ubound <- matrix(ubound, nrow, ncol, byrow)
			} else if (len == nrow * (ncol - 1) / 2) {
				if(byrow) {
					tmp <- matrix(as.numeric(NA), nrow, ncol)
					tmp[upper.tri(tmp)] <- ubound
					tmp[lower.tri(tmp)] <- t(tmp)[lower.tri(tmp)]
					ubound <- tmp
				} else {
					tmp <- matrix(as.numeric(NA), nrow, ncol)
					tmp[lower.tri(tmp)] <- ubound
					tmp[upper.tri(tmp)] <- t(tmp)[upper.tri(tmp)]
					ubound <- tmp
				}				
			} else {
				stop(paste(
					"Illegal number of elements (", len,
					") for ubound matrix in standardized matrix constructor", sep=""),
					call. = FALSE)
			}
		}
		retval <- callNextMethod(.Object, labels, values, free, lbound, ubound, name)
		return(retval)
	}
)

setMethod("omxVerifyMatrix", "StandMatrix",
	function(.Object) {
		callNextMethod(.Object)
		values <- .Object@values
		free <- .Object@free
		labels <- .Object@labels
		lbound <- .Object@lbound
		ubound <- .Object@ubound
		if (!all(values == t(values))) {
			stop(paste("Values matrix of standardized matrix", omxQuotes(.Object@name), 
				"is not symmetric!"), call. = FALSE)
		}
		if (!all(diag(values) == 1)) {
			stop(paste("Values matrix of standardized matrix", omxQuotes(.Object@name), 
				"is not 1's along the diagonal!"), call. = FALSE)
		}
		if (!all(free == t(free))) {
			stop(paste("Free matrix of standardized matrix", omxQuotes(.Object@name), 
				"is not symmetric!"), call. = FALSE)
		}
		if (!all(diag(free) == FALSE)) {
			stop(paste("Free matrix of standardized matrix", omxQuotes(.Object@name), 
				"is not fixed along the diagonal!"), call. = FALSE)
		}
		if (!all(labels == t(labels), na.rm = TRUE) && all(is.na(labels) == is.na(t(labels)))) {
			stop(paste("Labels matrix of standardized matrix", omxQuotes(.Object@name), 
				"is not symmetric!"), call. = FALSE)
		}
		if (!all(is.na(diag(labels)))) {
			stop(paste("Labels matrix of standardized matrix", omxQuotes(.Object@name), 
				"is not NA along the diagonal!"), call. = FALSE)
		}
		if (!all(lbound == t(lbound), na.rm = TRUE) && all(is.na(lbound) == is.na(t(lbound)))) {
			stop(paste("Lbound matrix of standardized matrix", omxQuotes(.Object@name), 
				"is not symmetric!"), call. = FALSE)
		}
		if (!all(is.na(diag(lbound)))) {
			stop(paste("Lbound matrix of standardized matrix", omxQuotes(.Object@name), 
				"is not NA along the diagonal!"), call. = FALSE)
		}
		if (!all(ubound == t(ubound), na.rm = TRUE) && all(is.na(ubound) == is.na(t(ubound)))) {
			stop(paste("Ubound matrix of standardized matrix", omxQuotes(.Object@name), 
				"is not symmetric!"), call. = FALSE)
		}
		if (!all(is.na(diag(ubound)))) {
			stop(paste("Ubound matrix of standardized matrix", omxQuotes(.Object@name), 
				"is not NA along the diagonal!"), call. = FALSE)
		}
	}
)
