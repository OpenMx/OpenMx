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


setClass(Class = "StandMatrix",
	representation = representation(),
	contains = "MxMatrix")

setMethod("imxSymmetricMatrix", "StandMatrix",
	function(.Object) { return(TRUE) }
)

setMethod("imxSquareMatrix", "StandMatrix",
	function(.Object) { return(TRUE) }
)

populateStandTriangle <- function(input, n, default, byrow, strname) {
	len <- length(input)
	if (len == 1) {
		output <- matrix(default, n, n)
		output[row(output) != col(output)] <- input
	} else if (len == n * n) {
		output <- matrix(input, n, n, byrow)
	} else if (len == n * (n - 1) / 2) {
		if(byrow) {
			output <- matrix(default, n, n)
			output[upper.tri(output)] <- input
			output[lower.tri(output)] <- t(output)[lower.tri(output)]
		} else {
			output <- matrix(default, n, n)
			output[lower.tri(output)] <- input
			output[upper.tri(output)] <- t(output)[upper.tri(output)]
		}				
	} else {
		stop(paste(
			"illegal number of elements (", len,
			") for '", strname, "' matrix in Standardized MxMatrix construction ", sep="",
			deparse(width.cutoff = 400L, imxLocateFunction("mxMatrix"))),
			call. = FALSE)
	}
	return(output)
}
	
setMethod("imxCreateMatrix", "StandMatrix",
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
			values <- populateStandTriangle(values, nrow, 1, byrow, 'values') 
		}
		if(condenseSlots && all.false(free) && all.na(labels)){
		  labels <- as.character(NA)
		  free <- FALSE
		}
    else{
  		if (is.vector(labels)) {
  			labels <- populateStandTriangle(labels, nrow, as.character(NA), byrow, 'labels')
  		}
  		if (is.vector(free)) {
  			free <- populateStandTriangle(free, nrow, FALSE, byrow, 'free')
  	}}
    if(condenseSlots && all.na(lbound)){lbound <- as.numeric(NA)}
		else{if (is.vector(lbound)) {
			lbound <- populateStandTriangle(lbound, nrow, as.numeric(NA), byrow, 'lbound')
		}}
		if(condenseSlots && all.na(ubound)){ubound <- as.numeric(NA)}
    else{if (is.vector(ubound)) {
			ubound <- populateStandTriangle(ubound, nrow, as.numeric(NA), byrow, 'ubound')
		}}
		return(callNextMethod(.Object, labels, values, free, lbound, ubound, nrow, ncol, byrow, name,
				      condenseSlots, joinKey, joinModel))
	}
)

setMethod("imxVerifyMatrix", "StandMatrix",
	function(.Object) {
		callNextMethod(.Object)
		values <- .Object@values
		free <- .Object@free
		labels <- .Object@labels
		lbound <- .Object@lbound
		ubound <- .Object@ubound
		if (!all(values == t(values))) {
			stop(paste("'values' matrix of Standardized MxMatrix", omxQuotes(.Object@name), 
				"is not symmetric in", 
				deparse(width.cutoff = 400L, imxLocateFunction("mxMatrix"))),
				call. = FALSE)
		}
		if (!all(diag(values) == 1)) {

		if (max(abs(diag(values) - 1)) < 1.0e-8) {
				stop(paste("'values' matrix of Standardized MxMatrix", omxQuotes(.Object@name), 
					"is very near, but not equal to, 1s along the diagonal in", 
					deparse(width.cutoff = 400L, imxLocateFunction("mxMatrix"))),
					call. = FALSE)
			} else {
				stop(paste("'values' matrix of Standardized MxMatrix", omxQuotes(.Object@name), 
					"is not 1s along the diagonal in", 
					deparse(width.cutoff = 400L, imxLocateFunction("mxMatrix"))),
					call. = FALSE)
			}
		}
		if (!all(free == t(free))) {
			stop(paste("'free' matrix of Standardized MxMatrix", omxQuotes(.Object@name), 
				"is not symmetric in", 
				deparse(width.cutoff = 400L, imxLocateFunction("mxMatrix"))),
				call. = FALSE)
		}
		if (!all(diag(free) == FALSE)) {
			stop(paste("'free' matrix of Standardized MxMatrix", omxQuotes(.Object@name), 
				"is not fixed along the diagonal in", 
				deparse(width.cutoff = 400L, imxLocateFunction("mxMatrix"))),
				call. = FALSE)
		}
		if (!all(labels == t(labels), na.rm = TRUE) && all(is.na(labels) == is.na(t(labels)))) {
			stop(paste("'labels' matrix of Standardized MxMatrix", omxQuotes(.Object@name), 
				"is not symmetric in", 
				deparse(width.cutoff = 400L, imxLocateFunction("mxMatrix"))),
				call. = FALSE)
		}
		if (!all(is.na(diag(labels)))) {
			stop(paste("'labels' matrix of Standardized MxMatrix", omxQuotes(.Object@name), 
				"is not NA along the diagonal in", 
				deparse(width.cutoff = 400L, imxLocateFunction("mxMatrix"))),
				call. = FALSE)
		}
		if (!all(lbound == t(lbound), na.rm = TRUE) && all(is.na(lbound) == is.na(t(lbound)))) {
			stop(paste("'lbound' matrix of Standardized MxMatrix", omxQuotes(.Object@name), 
				"is not symmetric in", 
				deparse(width.cutoff = 400L, imxLocateFunction("mxMatrix"))),
				call. = FALSE)
		}
		if (!all(is.na(diag(lbound)))) {
			stop(paste("'lbound' matrix of Standardized MxMatrix", omxQuotes(.Object@name), 
				"is not NA along the diagonal in", 
				deparse(width.cutoff = 400L, imxLocateFunction("mxMatrix"))),
				call. = FALSE)
		}
		if (!all(ubound == t(ubound), na.rm = TRUE) && all(is.na(ubound) == is.na(t(ubound)))) {
			stop(paste("'ubound' matrix of Standardized MxMatrix", omxQuotes(.Object@name), 
				"is not symmetric in", 
				deparse(width.cutoff = 400L, imxLocateFunction("mxMatrix"))),
				call. = FALSE)
		}
		if (!all(is.na(diag(ubound)))) {
			stop(paste("'ubound' matrix of Standardized MxMatrix", omxQuotes(.Object@name), 
				"is not NA along the diagonal in", 
				deparse(width.cutoff = 400L, imxLocateFunction("mxMatrix"))),
				call. = FALSE)
		}
	}
)
