#
#   Copyright 2007-20010 The OpenMx Project
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


setClass(Class = "SdiagMatrix",
	representation = representation(),
	contains = "MxMatrix")

setMethod("imxSymmetricMatrix", "SdiagMatrix",
	function(.Object) { return(FALSE) }
)

setMethod("imxSquareMatrix", "SdiagMatrix",
	function(.Object) { return(TRUE) }
)

populateSdiagTriangle <- function(input, n, default, byrow, strname) {
	len <- length(input)
	if (len == n * n) {
		output <- matrix(input, n, n, byrow)
	} else if (len == n * (n - 1) / 2 || len == 1) {
		if(byrow) {
			output <- matrix(default, n, n)
			output[upper.tri(output)] <- input
			output[lower.tri(output)] <- t(output)[lower.tri(output)]
			output[upper.tri(output, TRUE)] <- default
		} else {
			output <- matrix(default, n, n)
			output[lower.tri(output)] <- input
		}
	} else {
		stop(paste(
			"illegal number of elements (", len,
			") for '", strname, "' matrix of Subdiagonal MxMatrix construction", sep="",
			deparse(width.cutoff = 400L, imxLocateFunction("mxMatrix"))), call. = FALSE)
	}
	return(output)
}
	
setMethod("imxCreateMatrix", "SdiagMatrix",
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
			values <- populateSdiagTriangle(values, nrow, 0, byrow, 'values')
		}
		if(condenseSlots && all.false(free) && all.na(labels)){
		  labels <- as.character(NA)
		  free <- FALSE
		}
		else{
      if (is.vector(labels)) {
			  labels <- populateSdiagTriangle(labels, nrow, as.character(NA), byrow, 'labels')
		  }
		  if (is.vector(free)) {
			  free <- populateSdiagTriangle(free, nrow, FALSE, byrow, 'free')
		}}
    if(condenseSlots && all.na(lbound)){lbound <- as.numeric(NA)}
		else{if (is.vector(lbound)) {
			lbound <- populateSdiagTriangle(lbound, nrow, as.numeric(NA), byrow, 'lbound')
		}}
    if(condenseSlots && all.na(ubound)){ubound <- as.numeric(NA)}
		else{if (is.vector(ubound)) {
			ubound <- populateSdiagTriangle(ubound, nrow, as.numeric(NA), byrow, 'ubound')
		}}
		return(callNextMethod(.Object, labels, values, free, lbound, ubound, nrow, ncol, byrow, name,
				      condenseSlots, joinKey, joinModel))
	}
)

setMethod("imxVerifyMatrix", "SdiagMatrix",
	function(.Object) {
		callNextMethod(.Object)
		values <- .Object@values
		free <- .Object@free
		labels <- .Object@labels
		lbound <- .Object@lbound
		ubound <- .Object@ubound
		if (!all(values[upper.tri(values, TRUE)]  == 0)) {
			stop(paste("upper triangle or diagonal of 'values' matrix in Subdiagonal MxMatrix", omxQuotes(.Object@name), 
				"is not all zeros in", 
				deparse(width.cutoff = 400L, imxLocateFunction("mxMatrix"))),
				call. = FALSE)
		}
		if (any(free) && !all(free[upper.tri(free, TRUE)] == FALSE)) {
			stop(paste("upper triangle or diagonal of 'free' matrix in Subdiagonal MxMatrix", omxQuotes(.Object@name), 
				"is not all fixed in", 
				deparse(width.cutoff = 400L, imxLocateFunction("mxMatrix"))),
				call. = FALSE)
		}
		if (!all.na(labels) && !all(is.na(labels[upper.tri(labels, TRUE)]))) {
			stop(paste("upper triangle or diagonal of 'labels' matrix in Subdiagonal MxMatrix", omxQuotes(.Object@name), 
				"is not all NAs in", 
				deparse(width.cutoff = 400L, imxLocateFunction("mxMatrix"))),
				call. = FALSE)
		}
		if (!all.na(lbound) && !all(is.na(lbound[upper.tri(labels, TRUE)]))) {
			stop(paste("upper triangle or diagonal of 'lbound' matrix in Subdiagonal MxMatrix", omxQuotes(.Object@name), 
				"is not all NAs in", 
				deparse(width.cutoff = 400L, imxLocateFunction("mxMatrix"))),
				call. = FALSE)
		}
		if (!all.na(ubound) && !all(is.na(ubound[upper.tri(labels, TRUE)]))) {
			stop(paste("upper triangle or diagonal of 'ubound' matrix in Subdiagonal MxMatrix", omxQuotes(.Object@name), 
				"is not all NAs in", 
				deparse(width.cutoff = 400L, imxLocateFunction("mxMatrix"))),
				call. = FALSE)
		}
	}
)
