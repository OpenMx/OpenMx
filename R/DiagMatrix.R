#
#   Copyright 2007-2011 The OpenMx Project
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

setClass(Class = "DiagMatrix",
	representation = representation(),
	contains = "MxMatrix")

setMethod("imxSquareMatrix", "DiagMatrix",
	function(.Object) { return(TRUE) }
)

setMethod("imxCreateMatrix", "DiagMatrix",
	function(.Object, labels, values, free, lbound, ubound, nrow, ncol, byrow, name, ...) {
		if (nrow != ncol) {
			stop("Non-square matrix attempted for diagonal matrix constructor", call. = FALSE)
		}
		if (single.na(values)) {
			values <- 0
		}
		if (is.vector(values)) {
			tmp <- matrix(0, nrow, ncol)
			diag(tmp) <- values
			values <- tmp
		}
		if (is.vector(labels)) {
			tmp <- matrix(as.character(NA), nrow, ncol)
			diag(tmp) <- labels
			labels <- tmp
		}
		if (is.vector(free)) {
			tmp <- matrix(FALSE, nrow, ncol)
			diag(tmp) <- free
			free <- tmp
		}
		if (is.vector(lbound)) {
			tmp <- matrix(as.numeric(NA), nrow, ncol)
			diag(tmp) <- lbound
			lbound <- tmp
		}
		if (is.vector(ubound)) {
			tmp <- matrix(as.numeric(NA), nrow, ncol)
			diag(tmp) <- ubound
			ubound <- tmp
		}
		retval <- callNextMethod(.Object, labels, values, free, lbound, ubound, nrow, ncol, byrow, name, ...)
		return(retval)
	}
)

setMethod("imxVerifyMatrix", "DiagMatrix",
	function(.Object) {
		callNextMethod(.Object)
		values <- .Object@values
		free <- .Object@free
		if(nnzero(values[row(values) != col(values)]) > 0)
			{ stop(paste("Values matrix of", .Object@name, "is not a diagonal matrix.")) }
		if(any(free[row(free) != col(free)])) {
			{ stop(paste("Free matrix of", .Object@name, "has TRUE on non-diagonal.")) }
		}		
	}
)
