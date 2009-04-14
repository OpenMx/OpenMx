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


setClass(Class = "ZeroMatrix",
	representation = representation(),
	contains = "MxNonSymmetricMatrix")

setMethod("initialize", "ZeroMatrix",
	function(.Object, name, values, spec, nrow, ncol, byrow, free) {
		if (!single.na(values)) {
			warning("Ignoring values matrix for ZeroMatrix construction")
		}
		if (!single.na(spec)) {
			warning("Ignoring specification matrix for ZeroMatrix construction")
		}
		if (free) {
			warning("Ignoring \'free\' parameter for ZeroMatrix construction")
		}
		spec <- new("MxSparseMatrix", 0, nrow, ncol)
		values <- Matrix(0, nrow, ncol)
		return(callNextMethod(.Object, spec, values, name))
	}
)

setMethod("omxVerifyMatrix", "ZeroMatrix",
	function(.Object) {
		callNextMethod(.Object)	
		if(nnzero(.Object@spec) > 0) { 
			stop(paste("Specification matrix of zero matrix",
				omxQuotes(.Object@name), "is not empty."), call. = FALSE) 
		} 
		if(nnzero(.Object@values) > 0) { 
			stop(paste("Values matrix of zero matrix",
				omxQuotes(.Object@name), "is not empty"), call. = FALSE)
		} 
	}
)
