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
	function(.Object, name, values, spec, nrow, ncol, byrow, free) {
		if (nrow != ncol) {
			stop("Non-square matrix attempted for SymmMatrix constructor")
		}
		if (!single.na(values) && is.vector(values)) {
			if (length(values) == (nrow * (nrow + 1) / 2)) {
			    mvalues <- matrix(0, nrow, ncol)
			    mvalues[lower.tri(mvalues, diag = TRUE)] <- values
			    mvalues <- mvalues + t(mvalues) - diag(mvalues) * diag(nrow)
			    values <- mvalues
			} else if (length(values) == (nrow * ncol)) {
				values <- matrix(values, nrow, ncol)
			} else {
				stop("Invalid length of values matrix for SymmMatrix constructor")
			}
		}
		if (!single.na(spec) && is.vector(spec)) {
			if (length(spec) == (nrow * (nrow + 1) / 2)) {
			    mspec <- matrix(0, nrow, ncol)
			    mspec[lower.tri(mspec, diag = TRUE)] <- spec
			    mspec <- mspec + t(mspec) - diag(mspec) * diag(nrow)
			    spec <- mspec
			} else if (length(spec) == (nrow * ncol)) {
				spec <- matrix(spec, nrow, ncol)
			} else {
				stop("Invalid length of specification matrix for SymmMatrix constructor")
			}
		}
		if (is(spec, "MxSymmetricSparse")) {
		} else if (single.na(spec) && free) {
			spec <- new("MxSymmetricSparse", matrix(NA, nrow, ncol))
	    } else if (single.na(spec)){
			spec <- new("MxSymmetricSparse", 0, nrow, ncol)
		} else if (is(spec, "Matrix")) {
			spec <- new("MxSymmetricSparse", as.matrix(spec))
	    } else {
	    	spec <- new("MxSymmetricSparse", matrix(spec, nrow, ncol))
	    }
	    if (is(values, "MxSymmetricSparse")) {
	    } else if (single.na(values)) {
	    	values <- new("MxSymmetricSparse", matrix(0, nrow, ncol))
	    } else if (is(values, "Matrix")) {
			values <- new("MxSymmetricSparse", as.matrix(values))
		} else {
	    	values <- new("MxSymmetricSparse", matrix(values, nrow, ncol))
	    }
		retval <- callNextMethod(.Object, spec, values, name)
		return(retval)
	}
)

setMethod("omxVerifyMatrix", "SymmMatrix",
	function(.Object) {
		callNextMethod(.Object)
		verifySquare(.Object)
		values <- as.matrix(.Object@values)
		if (!all(values == t(values))) {
			stop(paste("Symmetric matrix",omxQuotes(.Object@name),"is not symmetric!"), 
				call.=FALSE)
		}
	}
)
