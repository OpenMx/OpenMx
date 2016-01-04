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


setClass(Class = "UnitMatrix",
	representation = representation(),
	contains = "MxMatrix")
	
setMethod("imxCreateMatrix", "UnitMatrix",
	  function(.Object, labels, values, free, lbound, ubound, nrow, ncol, byrow, name,
		    condenseSlots, joinKey, joinModel) {
		if (!single.na(values)) {
			warning("ignoring 'values' matrix for Unit MxMatrix construction in ",
			        deparse(width.cutoff = 400L, imxLocateFunction("mxMatrix")), 
                                call. = FALSE)
		}
		if (!single.na(labels)) {
			warning("ignoring 'labels' matrix for Unit MxMatrix construction in ",
			        deparse(width.cutoff = 400L, imxLocateFunction("mxMatrix")), 
                                call. = FALSE)
		}
		if (!(length(free) == 1 && free == FALSE)) {
			warning("ignoring 'free' matrix for Unit MxMatrix construction in ",
			        deparse(width.cutoff = 400L, imxLocateFunction("mxMatrix")), 
                                call. = FALSE)
		}
		if (!single.na(lbound)) {
			warning("ignoring 'lbound' matrix for Unit MxMatrix construction in ",
			        deparse(width.cutoff = 400L, imxLocateFunction("mxMatrix")), 
                                 call. = FALSE)
		}
		if (!single.na(ubound)) {
			warning("ignoring 'ubound' matrix for Unit MxMatrix construction in ", 
			        deparse(width.cutoff = 400L, imxLocateFunction("mxMatrix")),
                                call. = FALSE)
		}
		labels <- matrix(as.character(NA), ifelse(condenseSlots,1,nrow), ifelse(condenseSlots,1,ncol))
		values <- matrix(1, nrow, ncol)
		free <- matrix(FALSE, ifelse(condenseSlots,1,nrow), ifelse(condenseSlots,1,ncol))
		lbound <- matrix(as.numeric(NA), ifelse(condenseSlots,1,nrow), ifelse(condenseSlots,1,ncol))
		ubound <- matrix(as.numeric(NA), ifelse(condenseSlots,1,nrow), ifelse(condenseSlots,1,ncol))
		return(callNextMethod(.Object, labels, values, free, lbound, ubound, nrow, ncol, byrow, name,
				      condenseSlots, joinKey, joinModel))
	}
)

setMethod("imxVerifyMatrix", "UnitMatrix",
	function(.Object) {
		callNextMethod(.Object)		
		if(!all(.Object@free == FALSE)) { 
			stop(paste("'free' matrix of Unit MxMatrix", 
				omxQuotes(.Object@name), "has a free parameter in "), 
				deparse(width.cutoff = 400L, imxLocateFunction("mxMatrix")),
				call.=FALSE)
		} 
		if(nnzero(.Object@values - 1) > 0) { 
			stop(paste("'values' matrix of Unit MxMatrix",
				omxQuotes(.Object@name), "has non unit entries in "), 
				deparse(width.cutoff = 400L, imxLocateFunction("mxMatrix")),
				call.=FALSE)
		} 
	}
)
