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


setClass(Class = "SymmMatrix",
	representation = representation(),
	contains = "MxMatrix")

setMethod("imxSymmetricMatrix", "SymmMatrix",
	function(.Object) { return(TRUE) }
)

setMethod("imxSquareMatrix", "SymmMatrix",
	function(.Object) { return(TRUE) }
)

populateSymmTriangle <- function(input, n, default, byrow, strname) {	
	len <- length(input)
	if (len == n * n || len == 1) {
		output <- matrix(input, n, n, byrow)
	} else if (len == n * (n + 1) / 2) {
		if(byrow) {
			output <- matrix(default, n, n)
			output[upper.tri(output, TRUE)] <- input
			output[lower.tri(output)] <- t(output)[lower.tri(output)]
		} else {
			output <- matrix(default, n, n)
			output[lower.tri(output, TRUE)] <- input
			output[upper.tri(output)] <- t(output)[upper.tri(output)]
		}				
	} else {
		stop(paste(
			"illegal number of elements (", len,
			") for '", strname, "' matrix in Symmetric MxMatrix construction", sep=""),
			deparse(width.cutoff = 400L, imxLocateFunction("mxMatrix")),
			call. = FALSE)
	}
	return(output)
}
	
	
setMethod("imxCreateMatrix", "SymmMatrix",
	  function(.Object, labels, values, free, lbound, ubound, nrow, ncol, byrow, name, condenseSlots,
		   joinKey, joinModel) {
		if (nrow != ncol) {
			stop(paste("non-square MxMatrix attempted in 'nrow' and 'ncol' arguments to",
			     deparse(width.cutoff = 400L, imxLocateFunction("mxMatrix"))), 
                             call. = FALSE)
		}
		if (single.na(values)) {
			values <- 0
		}
		if (is.vector(values)) {
			values <- populateSymmTriangle(values, nrow, 0, byrow, 'values')
		}
		if(condenseSlots && all.false(free) && all.na(labels)){
		  labels <- as.character(NA)
		  free <- FALSE
		}
    else{
  		if (is.vector(labels)) {
  			labels <- populateSymmTriangle(labels, nrow, as.character(NA), byrow, 'labels')
  		}
  		if (is.vector(free)) {
  			free <- populateSymmTriangle(free, nrow, FALSE, byrow, 'free')
		}}
    if(condenseSlots && all.na(lbound)){lbound <- as.numeric(NA)}
		else{if (is.vector(lbound)) {
			lbound <- populateSymmTriangle(lbound, nrow, as.numeric(NA), byrow, 'lbound')
		}}
		if(condenseSlots && all.na(ubound)){ubound <- as.numeric(NA)}
    else{if (is.vector(ubound)) {
			ubound <- populateSymmTriangle(ubound, nrow, as.numeric(NA), byrow, 'ubound')
		}}
		return(callNextMethod(.Object, labels, values, free, lbound, ubound, nrow, ncol, byrow, name,
				      condenseSlots, joinKey, joinModel))
	}
)

setMethod("imxVerifyMatrix", "SymmMatrix",
	function(.Object) {
		callNextMethod(.Object)
		values <- .Object@values
		free <- .Object@free
		labels <- .Object@labels
		lbound <- .Object@lbound
		ubound <- .Object@ubound
		mask   <- is.na(values[upper.tri(values)])
		if (any(is.na(t(values)[upper.tri(values)]) != mask)) {
			stop(paste("NAs in 'values' matrix of Symmetric MxMatrix", omxQuotes(.Object@name), 
				"are not symmetric in "), 
				deparse(width.cutoff = 400L, imxLocateFunction("mxMatrix")),
				call. = FALSE)
		}
		if (any(values[upper.tri(values)][!mask]  != t(values)[upper.tri(values)][!mask])) {
			stop(paste("'values' matrix of Symmetric MxMatrix", omxQuotes(.Object@name), 
				"is not symmetric in "), 				deparse(width.cutoff = 400L, imxLocateFunction("mxMatrix")),
				call. = FALSE)
		}
		if (!all(free == t(free))) {
			stop(paste("'free' matrix of Symmetric MxMatrix", omxQuotes(.Object@name), 
				"is not symmetric in "), 
				deparse(width.cutoff = 400L, imxLocateFunction("mxMatrix")),
				call. = FALSE)
		}
		if (!all(labels == t(labels), na.rm = TRUE) && all(is.na(labels) == is.na(t(labels)))) {
			stop(paste("'labels' matrix of Symmetric MxMatrix", omxQuotes(.Object@name), 
				"is not symmetric in "), 
				deparse(width.cutoff = 400L, imxLocateFunction("mxMatrix")),
				call. = FALSE)
		}
		if (!all(lbound == t(lbound), na.rm = TRUE) && all(is.na(lbound) == is.na(t(lbound)))) {
			stop(paste("'lbound' matrix of Symmetric MxMatrix", omxQuotes(.Object@name), 
				"is not symmetric in "), 
				deparse(width.cutoff = 400L, imxLocateFunction("mxMatrix")),
				call. = FALSE)
		}
		if (!all(ubound == t(ubound), na.rm = TRUE) && all(is.na(ubound) == is.na(t(ubound)))) {
			stop(paste("'ubound' matrix of Symmetric MxMatrix", omxQuotes(.Object@name), 
				"is not symmetric in "), 
				deparse(width.cutoff = 400L, imxLocateFunction("mxMatrix")),
				call. = FALSE)
		}
	}
)
