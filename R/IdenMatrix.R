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


setClass(Class = "IdenMatrix",
	representation = representation(),
	contains = "MxNonSymmetricMatrix")
	
setMethod("initialize", "IdenMatrix",
	function(.Object, name, values, free, labels, nrow, ncol, byrow) {
		if (!single.na(values)) {
			warning("Ignoring values matrix for IdenMatrix construction", call. = FALSE)
		}
		if (!single.na(labels)) {
			warning("Ignoring labels matrix for IdenMatrix construction", call. = FALSE)
		}
		if (!(length(free) == 1 && free == FALSE)) {
			warning("Ignoring free matrix for IdenMatrix construction", call. = FALSE)
		}		
		if (nrow != ncol) {
			stop("Non-square matrix attempted for IdenMatrix constructor", call. = FALSE)
		}
		labels <- matrix("", nrow, ncol)
		values <- matrix(diag(nrow = nrow), nrow, ncol)
		free <- matrix(FALSE, nrow, ncol)
		return(callNextMethod(.Object, labels, values, free, name))
	}
)

setMethod("omxVerifyMatrix", "IdenMatrix",
	function(.Object) {
		callNextMethod(.Object)
		verifySquare(.Object)
		if(!all(.Object@free == FALSE)) {
			stop(paste("Free matrix of iden matrix", 
				omxQuotes(.Object@name), "has a free parameter"), call.=FALSE)
		} 
		if(!suppressWarnings(all(.Object@values == diag(nrow(.Object@values))))) {
			stop(paste("Values matrix of iden matrix",
				omxQuotes(.Object@name), "is not the identity matrix"), call.=FALSE)
		}
	}
)
