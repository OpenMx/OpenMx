#
#   Copyright 2007-2009 The OpenMx Project
#
#   Licensed under the Apache License, Version 2.0 (the "License");
#   you may not use this file except in compliance with the License.
#   You may obtain a copy of the License at
# 
#        http://www.apache.org/licenses/LICENSE-2.0
# 
#    Unless required by applicable law or agreed to in writing, software
#    distributed under the License is distributed on an "AS IS" BASIS,
#    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#   See the License for the specific language governing permissions and
#   limitations under the License.


setClass(Class = "IdenMatrix",
	representation = representation(),
	contains = "MxNonSymmetricMatrix")
	
setMethod("initialize", "IdenMatrix",
	function(.Object, name, values, spec, nrow, ncol, byrow, free) {
		if (!single.na(values)) {
			warning("Ignoring values matrix for IdenMatrix construction")
		}
		if (!single.na(spec)) {
			warning("Ignoring specification matrix for IdenMatrix construction")
		}
		if (free) {
			warning("Ignoring \'free\' parameter for IdenMatrix construction")
		}		
		if (nrow != ncol) {
			stop("Non-square matrix attempted for IdenMatrix constructor")
		}
		spec <- new("MxSparseMatrix", 0, nrow, ncol)
		values <- Matrix(diag(nrow))
		return(callNextMethod(.Object, spec, values, name))
	}
)

setMethod("omxVerifyMatrix", "IdenMatrix",
	function(.Object) {
		callNextMethod(.Object)
		verifySquare(.Object)
		if(nnzero(.Object@spec) > 0) { stop("Specification matrix is not empty") } 
		if(!suppressWarnings(all(.Object@values == diag(nrow(.Object@values))))) 
			{ stop("Values matrix is not the identity matrix") }
	}
)
