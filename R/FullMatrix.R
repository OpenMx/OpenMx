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


setClass(Class = "FullMatrix",
	representation = representation(),
	contains = "MxNonSymmetricMatrix")

setMethod("initialize", "FullMatrix",
	function(.Object, name, values, specification, nrow, ncol, byrow, free) {
		if (is(specification, "MxSparseMatrix")) {
		} else if (single.na(specification) && free) {
			specification <- new("MxSparseMatrix", matrix(NA, nrow, ncol))
	    } else if (single.na(specification)) {
			specification <- new("MxSparseMatrix", 0, nrow, ncol)
	    } else if (is(specification, "Matrix")) {
	    	specification <- new("MxSparseMatrix", as.matrix(specification))
	    } else if (is.matrix(specification)) {
	    	specification <- new("MxSparseMatrix", specification)
	    } else if(is.vector(specification)) {
	    	specification <- new("MxSparseMatrix", matrix(specification, nrow, ncol))
	    } else {
	    	specification <- new("MxSparseMatrix", specification, nrow, ncol)
	    }
		if (is(values, "Matrix")) {
		} else if (is.matrix(values)) {
			values <- Matrix(values)
	    } else if (single.na(values)) {
	    	values <- Matrix(0, nrow, ncol)
	    } else {
	    	values <- Matrix(values, nrow, ncol)
	    } 
		retval <- callNextMethod(.Object, specification, values, name) 
		return(retval)
	}
)

setMethod("omxVerifyMatrix", "FullMatrix",
	function(.Object) {
		callNextMethod(.Object)
	}
)
