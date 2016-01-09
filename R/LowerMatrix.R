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


setClass(Class = "LowerMatrix",
	representation = representation(),
	contains = "MxMatrix")

setMethod("imxSymmetricMatrix", "LowerMatrix",
	function(.Object) { return(FALSE) }
)

setMethod("imxSquareMatrix", "LowerMatrix",
	function(.Object) { return(TRUE) }
)

populateLowerTriangle <- function(input, n, default, byrow, strname) {
	len <- length(input)
	if (len == n * n) {
		output <- matrix(input, n, n, byrow)
	} else if (len == n * (n + 1) / 2 || len == 1) {
		if(byrow) {
			output <- matrix(default, n, n)
			output[upper.tri(output, TRUE)] <- input
			output[lower.tri(output)] <- t(output)[lower.tri(output)]
			output[upper.tri(output)] <- default
		} else {
			output <- matrix(default, n, n)
			output[lower.tri(output, TRUE)] <- input
		}			
	} else {
		stop(paste(
			"illegal number of elements (", len,
			") for '", strname, "' matrix in Lower MxMatrix construction ", sep="",
			deparse(width.cutoff = 400L, imxLocateFunction("mxMatrix"))), call. = FALSE)
	}
	return(output)
}
	
setMethod("imxCreateMatrix", "LowerMatrix",
	  function(.Object, labels, values, free, lbound, ubound, nrow, ncol, byrow, name,
		    condenseSlots, joinKey, joinModel) {
		if (nrow != ncol) {
			stop(paste("non-square MxMatrix attempted in 'nrow' and 'ncol' arguments to",
			     deparse(width.cutoff = 400L, imxLocateFunction("mxMatrix"))), 
                             call. = FALSE)
		}
		if (single.na(values)) {
			values <- 0
		}
		if (is.vector(values)) {
			values <- populateLowerTriangle(values, nrow, 0, byrow, 'values')
		}
		if(condenseSlots && all.false(free) && all.na(labels)){
		  labels <- as.character(NA)
		  free <- FALSE
		}
    else{
  		if (is.vector(labels)) {
  			labels <- populateLowerTriangle(labels, nrow, as.character(NA), byrow, 'labels')
  		}
  		if (is.vector(free)) {
  			free <- populateLowerTriangle(free, nrow, FALSE, byrow, 'free')
  	}}
    if(condenseSlots && all.na(lbound)){lbound <- as.numeric(NA)}
		else{if (is.vector(lbound)) {
			lbound <- populateLowerTriangle(lbound, nrow, as.numeric(NA), byrow, 'lbound')
		}}
		if(condenseSlots && all.na(ubound)){ubound <- as.numeric(NA)}
		else{if (is.vector(ubound)) {
			ubound <- populateLowerTriangle(ubound, nrow, as.numeric(NA), byrow, 'ubound')
		}}
		return(callNextMethod(.Object, labels, values, free, lbound, ubound, nrow, ncol, byrow, name,
				      condenseSlots, joinKey, joinModel))
	}
)

setMethod("imxVerifyMatrix", "LowerMatrix",
	function(.Object) {
		callNextMethod(.Object)
    #Do all of these slots really need to be copied?:
		values <- .Object@values
		free <- .Object@free
		labels <- .Object@labels
		lbound <- .Object@lbound
		ubound <- .Object@ubound
		if (!all(values[upper.tri(values)]  == 0)) {
			stop(paste("upper triangle of 'values' matrix in Lower MxMatrix", omxQuotes(.Object@name), 
				"is not all zeros in", 
				deparse(width.cutoff = 400L, imxLocateFunction("mxMatrix"))),
				call. = FALSE)
		}
		if (any(free) && !all(free[upper.tri(free)] == FALSE)) {
			stop(paste("upper triangle of 'free' matrix in Lower MxMatrix", omxQuotes(.Object@name), 
				"is not all fixed in", 
				deparse(width.cutoff = 400L, imxLocateFunction("mxMatrix"))),
				call. = FALSE)
		}
		if (!all.na(labels) && !all(is.na(labels[upper.tri(labels)]))) {
			stop(paste("upper triangle of 'labels' matrix in Lower MxMatrix", omxQuotes(.Object@name), 
				"is not all NAs in", 
				deparse(width.cutoff = 400L, imxLocateFunction("mxMatrix"))),
				call. = FALSE)
		}
		if (!all.na(lbound) && !all(is.na(lbound[upper.tri(lbound)]))) {
			stop(paste("upper triangle of 'lbound' matrix in Lower MxMatrix", omxQuotes(.Object@name), 
				"is not all NAs in", 
				deparse(width.cutoff = 400L, imxLocateFunction("mxMatrix"))),
				call. = FALSE)
		}
		if (!all.na(ubound) && !all(is.na(ubound[upper.tri(ubound)]))) {
			stop(paste("upper triangle of 'ubound' matrix in Lower MxMatrix", omxQuotes(.Object@name), 
				"is not all NAs in", 
				deparse(width.cutoff = 400L, imxLocateFunction("mxMatrix"))),
				call. = FALSE)
		}
	}
)
