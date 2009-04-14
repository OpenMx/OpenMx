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
	function(.Object, name, values, spec, nrow, ncol, byrow, free) {
		if (nrow != ncol) {
			stop("Non-square matrix attempted for DiagMatrix constructor")
		}
	    if (single.na(spec) && free) {
	    	spec <- diag(nrow)
	    	spec[spec == 1] <- NA
			spec <- new("MxSparseMatrix", spec)
		} else if (single.na(spec)) {
			spec <- new("MxSparseMatrix", 0, nrow, ncol)
 		} else if (is(spec, "Matrix")) {
			spec <- new("MxSparseMatrix", as.matrix(spec))
		} else if (is.matrix(spec)) {
			spec <- new("MxSparseMatrix", spec)
		} else if (is.vector(spec)) {
			tmp <- matrix(0, nrow, nrow)
			diag(tmp) <- spec
			spec <- new("MxSparseMatrix", tmp)
		} else {
			spec <- new("MxSparseMatrix", spec)
		}
		if (single.na(values)) {
			values <- Matrix(0, nrow, ncol)
		} else if (is.matrix(values)) { 
			values <- Matrix(values)
		} else if (is.vector(values)) {
			values <- Matrix(values * diag(nrow))
		} else {
			values <- Matrix(values)
		}
		retval <- callNextMethod(.Object, spec, values, name) 
		return(retval)
	}
)

setMethod("omxVerifyMatrix", "DiagMatrix",
	function(.Object) {
		callNextMethod(.Object)
		verifySquare(.Object)
		values <- .Object@values
		if(suppressWarnings(nnzero(values - diag(diag(values,
			nrow = nrow(values), ncol = ncol(values))))) > 0)
			{ stop(paste("Values matrix of", .Object@name, "is not a diagonal matrix")) }
	}
)
