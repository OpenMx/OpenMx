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

setClass(Class = "DiagMatrix",
	representation = representation(),
	contains = "MxMatrix")

setMethod("imxSquareMatrix", "DiagMatrix",
	function(.Object) { return(TRUE) }
)

setMethod("imxCreateMatrix", "DiagMatrix",
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
			tmp <- matrix(0, nrow, ncol)
			diag(tmp) <- values
			values <- tmp
		}
		if(condenseSlots && all.na(labels) && all.false(free)){
      labels <- as.character(NA)
		  free <- FALSE
		}
    else{
      if(is.vector(labels)) {
		    tmp <- matrix(as.character(NA), nrow, ncol)
		    diag(tmp) <- labels
		    labels <- tmp
      }
	  if(is.vector(free)) {
		  tmp <- matrix(FALSE, nrow, ncol)
		  diag(tmp) <- free
		  free <- tmp
		}}
    if(condenseSlots && all.na(lbound)){lbound <- as.numeric(NA)}
		else{if(is.vector(lbound)) {
		  tmp <- matrix(as.numeric(NA), nrow, ncol)
		  diag(tmp) <- lbound
		  lbound <- tmp
		}}
    if(condenseSlots && all.na(ubound)){ubound <- as.numeric(NA)}
		else{if(is.vector(ubound)) {
		  tmp <- matrix(as.numeric(NA), nrow, ncol)
		  diag(tmp) <- ubound
		  ubound <- tmp
		}}
    if(exists("tmp")){rm(tmp)}
		return(callNextMethod(.Object, labels, values, free, lbound, ubound, nrow, ncol, byrow, name,
				      condenseSlots, joinKey, joinModel))
	}
)

setMethod("imxVerifyMatrix", "DiagMatrix",
	function(.Object) {
		callNextMethod(.Object)
		values <- .Object@values
		if(nnzero(values[row(values) != col(values)]) > 0)
			{ stop(paste("'values' matrix of", .Object@name, "is not a diagonal matrix in",
			deparse(width.cutoff = 400L, imxLocateFunction("mxMatrix"))),
			call. = FALSE) }
    rm(values)
		free <- .Object@free
		if(any(free) && any(free[row(free) != col(free)])) {
			{ stop(paste("'free' matrix of", .Object@name, "has TRUE on non-diagonal in",
			deparse(width.cutoff = 400L, imxLocateFunction("mxMatrix"))),
			call. = FALSE) }
		}		
	}
)
