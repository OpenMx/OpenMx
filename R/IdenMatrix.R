#
#   Copyright 2007-2016 The OpenMx Project
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
	  function(.Object, labels, values, free, lbound, ubound, nrow, ncol, byrow, name,
		    condenseSlots, joinKey, joinModel) {
		if (!single.na(values) && !is.null(values)) {
			warning("ignoring 'values' matrix for Identity MxMatrix construction in ",
			        deparse(width.cutoff = 400L, imxLocateFunction("mxMatrix")), 
                                call. = FALSE)
		}
		if (!single.na(labels) && !is.null(labels)) {
			warning("ignoring 'labels' matrix for Identity MxMatrix construction in ",
			        deparse(width.cutoff = 400L, imxLocateFunction("mxMatrix")), 
                                call. = FALSE)
		}
		if (!(length(free) == 1 && free == FALSE) && !is.null(free)) {
			warning("ignoring 'free' matrix for Identity MxMatrix construction in ",
			        deparse(width.cutoff = 400L, imxLocateFunction("mxMatrix")),
                                 call. = FALSE)
		}
		if (!single.na(lbound) && !is.null(lbound)) {
			warning("ignoring 'lbound' matrix for Identity MxMatrix construction in ",
			        deparse(width.cutoff = 400L, imxLocateFunction("mxMatrix")), 
                                call. = FALSE)
		}
		if (!single.na(ubound) && !is.null(ubound)) {
			warning("ignoring 'ubound' matrix for Identity MxMatrix construction in ",
			        deparse(width.cutoff = 400L, imxLocateFunction("mxMatrix")), 
                                call. = FALSE)
		}
		if (nrow != ncol) {
			stop("non-square MxMatrix attempted for Identity MxMatrix construction in ",
			     deparse(width.cutoff = 400L, imxLocateFunction("mxMatrix")), 
                             call. = FALSE)
		}
		values <- matrix(diag(nrow = nrow), nrow, ncol)
		labels <- NULL
		free <- NULL
		lbound <- NULL
		ubound <- NULL
		if (!condenseSlots) {
			labels <- matrix(as.character(NA), nrow, ncol)
			free <- matrix(FALSE, nrow, ncol)
			lbound <- matrix(as.numeric(NA), nrow, ncol)
			ubound <- matrix(as.numeric(NA), nrow, ncol)
		}
		return(callNextMethod(.Object, labels, values, free, lbound, ubound, nrow, ncol, byrow, name,
				      condenseSlots, joinKey, joinModel))
	}
)

setMethod("imxVerifyMatrix", "IdenMatrix",
	function(.Object) {
		callNextMethod(.Object)
		if (!is.null(.Object@free) && !all(.Object@free == FALSE)) {
			stop(paste("'free' matrix of Identity MxMatrix", 
				omxQuotes(.Object@name), "has a free parameter in", 
				deparse(width.cutoff = 400L, imxLocateFunction("mxMatrix"))),
				call.=FALSE)
		}
		if(!suppressWarnings(all(.Object@values == diag(nrow(.Object@values))))) {
			stop(paste("'values' matrix of Identity MxMatrix",
				omxQuotes(.Object@name), "is not the identity matrix in", 
				deparse(width.cutoff = 400L, imxLocateFunction("mxMatrix"))),
				call.=FALSE)
		}
	}
)

