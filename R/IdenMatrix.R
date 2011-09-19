#
#   Copyright 2007-2010 The OpenMx Project
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
	contains = "MxMatrix")

setMethod("imxSquareMatrix", "IdenMatrix",
	function(.Object) { return(TRUE) }
)	

setMethod("imxCreateMatrix", "IdenMatrix",
	function(.Object, labels, values, free, lbound, ubound, nrow, ncol, byrow, name, ...) {
		if (!single.na(values)) {
			warning("Ignoring values matrix for identity matrix construction ",
			        deparse(width.cutoff = 400L, imxLocateFunction("mxMatrix")), 
                                call. = FALSE)
		}
		if (!single.na(labels)) {
			warning("Ignoring labels matrix for identity matrix construction ",
			        deparse(width.cutoff = 400L, imxLocateFunction("mxMatrix")), 
                                call. = FALSE)
		}
		if (!(length(free) == 1 && free == FALSE)) {
			warning("Ignoring free matrix for identity matrix construction ",
			        deparse(width.cutoff = 400L, imxLocateFunction("mxMatrix")),
                                 call. = FALSE)
		}
		if (!single.na(lbound)) {
			warning("Ignoring lbound matrix for identity matrix construction ",
			        deparse(width.cutoff = 400L, imxLocateFunction("mxMatrix")), 
                                call. = FALSE)
		}
		if (!single.na(ubound)) {
			warning("Ignoring ubound matrix for identity matrix construction ",
			        deparse(width.cutoff = 400L, imxLocateFunction("mxMatrix")), 
                                call. = FALSE)
		}
		if (nrow != ncol) {
			stop("Non-square matrix attempted for identity matrix constructor ",
			     deparse(width.cutoff = 400L, imxLocateFunction("mxMatrix")), 
                             call. = FALSE)
		}
		labels <- matrix(as.character(NA), nrow, ncol)
		values <- matrix(diag(nrow = nrow), nrow, ncol)
		free <- matrix(FALSE, nrow, ncol)
		lbound <- matrix(as.numeric(NA), nrow, ncol)
		ubound <- matrix(as.numeric(NA), nrow, ncol)
		return(callNextMethod(.Object, labels, values, free, lbound, ubound, nrow, ncol, byrow, name, ...))
	}
)

setMethod("imxVerifyMatrix", "IdenMatrix",
	function(.Object) {
		callNextMethod(.Object)
		if(!all(.Object@free == FALSE)) {
			stop(paste("Free matrix of identity matrix", 
				omxQuotes(.Object@name), "has a free parameter", 
				deparse(width.cutoff = 400L, imxLocateFunction("mxMatrix"))),
				call.=FALSE)
		}
		if(!suppressWarnings(all(.Object@values == diag(nrow(.Object@values))))) {
			stop(paste("Values matrix of identity matrix",
				omxQuotes(.Object@name), "is not the identity matrix", 
				deparse(width.cutoff = 400L, imxLocateFunction("mxMatrix"))),
				call.=FALSE)
		}
	}
)
