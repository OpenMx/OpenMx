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
	function(.Object, name, values, specification, nrow, ncol, byrow, free) {
		if (nrow != ncol) {
			stop("Non-square matrix attempted for DiagMatrix constructor")
		}
	    if (single.na(specification) && free) {
	    	specification <- diag(nrow)
	    	specification[specification == 1] <- NA
			specification <- new("MxSparseMatrix", specification)
		} else if (single.na(specification)) {
			specification <- new("MxSparseMatrix", 0, nrow, ncol)
 		} else if (is(specification, "Matrix")) {
			specification <- new("MxSparseMatrix", as.matrix(specification))
		} else if (is.matrix(specification)) {
			specification <- new("MxSparseMatrix", specification)
		} else if (is.vector(specification)) {
			tmp <- matrix(0, nrow, nrow)
			diag(tmp) <- specification
			specification <- new("MxSparseMatrix", tmp)
		} else {
			specification <- new("MxSparseMatrix", specification)
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
		retval <- callNextMethod(.Object, specification, values, name) 
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
