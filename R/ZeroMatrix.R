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
	contains = "MxMatrix")

setMethod("initialize", "ZeroMatrix",
	function(.Object, name, values, free, labels, lbound, ubound, nrow, ncol, byrow) {
		if (!single.na(values)) {
			warning("Ignoring values matrix for zero matrix constructor", call. = FALSE)
		}
		if (!single.na(labels)) {
			warning("Ignoring labels matrix for zero matrix constructor", call. = FALSE)
		}
		if (!(length(free) == 1 && free == FALSE)) {
			warning("Ignoring free matrix for zero matrix constructor", call. = FALSE)
		}
		if (!single.na(lbound)) {
			warning("Ignoring lbound matrix for zero matrix constructor", call. = FALSE)
		}
		if (!single.na(ubound)) {
			warning("Ignoring ubound matrix for zero matrix constructor", call. = FALSE)
		}
		labels <- matrix(as.character(NA), nrow, ncol)
		values <- matrix(0, nrow, ncol)
		free <- matrix(FALSE, nrow, ncol)
		lbound <- matrix(as.numeric(NA), nrow, ncol)
		ubound <- matrix(as.numeric(NA), nrow, ncol)
		return(callNextMethod(.Object, labels, values, free, lbound, ubound, name))
	}
)

setMethod("omxVerifyMatrix", "ZeroMatrix",
	function(.Object) {
		callNextMethod(.Object)	
		if(!all(.Object@free == FALSE)) { 
			stop(paste("Free matrix of zero matrix", 
				omxQuotes(.Object@name), "has a free parameter"), call.=FALSE)
		} 
		if(nnzero(.Object@values) > 0) { 
			stop(paste("Values matrix of zero matrix",
				omxQuotes(.Object@name), "has non zero entries"), call.=FALSE)
		} 
	}
)
