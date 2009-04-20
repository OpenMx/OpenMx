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

setClass(Class = "DiagMatrix",
	representation = representation(),
	contains = "MxNonSymmetricMatrix")

setMethod("initialize", "DiagMatrix",
	function(.Object, name, values, free, labels, nrow, ncol, byrow) {
		if (nrow != ncol) {
			stop("Non-square matrix attempted for DiagMatrix constructor", call. = FALSE)
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
			tmp <- matrix("", nrow, ncol)
			diag(tmp) <- labels
			labels <- tmp
		}
		if (is.vector(free)) {
			tmp <- matrix(FALSE, nrow, ncol)
			diag(tmp) <- free
			free <- tmp
		}
		retval <- callNextMethod(.Object, labels, values, free, name)
		return(retval)
	}
)

setMethod("omxVerifyMatrix", "DiagMatrix",
	function(.Object) {
		callNextMethod(.Object)
		verifySquare(.Object)
		values <- .Object@values
		free <- .Object@free
		if(suppressWarnings(nnzero(values - diag(diag(values,
			nrow = nrow(values), ncol = ncol(values))))) > 0)
			{ stop(paste("Values matrix of", .Object@name, "is not a diagonal matrix.")) }
		if(any(free[row(free) != col(free)])) {
			{ stop(paste("Free matrix of", .Object@name, "has TRUE on non-diagonal.")) }
		}		
	}
)
