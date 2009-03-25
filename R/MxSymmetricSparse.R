#
#   Copyright 2007-2009 The OpenMx Project
#
#   Licensed under the Apache License, Version 2.0 (the "License");
#   you may not use this file except in compliance with the License.
#   You may obtain a copy of the License at
# 
#        http://www.apache.org/licenses/LICENSE-2.0
# 
#    Unless required by applicable law or agreed to in writing, software
#    distributed under the License is distributed on an "AS IS" BASIS,
#    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#   See the License for the specific language governing permissions and
#   limitations under the License.


setClass(Class = "MxSymmetricSparse",
	representation = representation(
		rowVector = "vector",
		colVector = "vector",
		dataVector = "vector",
		nrow = "numeric",
		ncol = "numeric",
		dimNames = "list"),
	prototype = prototype(
		rowVector = vector(mode="numeric",0),
		colVector = vector(mode="numeric",0),
		dataVector = vector(mode="character",0),
		dimNames = list(),
		nrow = 1, ncol = 1))

setMethod("initialize", "MxSymmetricSparse",
	function(.Object, data = NA, nrow = 1, ncol = 1, dimnames = NULL) {
		if (is.matrix(data)) {
		    test <- all(data == t(data) | is.na(data))
		    if (is.na(test) || !test) {
				stop("Input matrix is not symmetric")
		    }
		    data[lower.tri(data)] <- 0
			tdata <- t(data)
			nonZero <- (data != 0) | is.na(data)
			rowMatrix <- t(row(data) * nonZero)
			colMatrix <- t(col(data) * nonZero)
			rowVector <- as.vector(rowMatrix)
			colVector <- as.vector(colMatrix)
			.Object@dataVector <- tdata[tdata != 0]
			.Object@rowVector <- rowVector[rowVector > 0]
			.Object@colVector <- colVector[colVector > 0]
			.Object@nrow <- nrow(data)
			.Object@ncol <- ncol(data)
		} else {
			.Object@nrow <- nrow
			.Object@ncol <- ncol
	    }
	    if (!is.null(dimnames)) {
	    	.Object@dimNames <- dimnames
	    }
		return(.Object)
	}
)

setMethod("print", "MxSymmetricSparse", function(x,...) { 
  args <- list(...)
  use.quotes <- args[['use.quotes']]
  if (is.null(use.quotes)) {
    use.quotes <- FALSE
  }
  omxDisplayMxSymmetricSparse(x, use.quotes) 
})

setMethod("show", "MxSymmetricSparse", function(object) { 
  omxDisplayMxSymmetricSparse(object, FALSE) 
})

omxDisplayMxSymmetricSparse <- function(mxMatrix, use.quotes) {
   matrix <- as.matrix(mxMatrix)
   if (use.quotes) {
      matrix <- apply(matrix, c(1,2), function(x) {
    if(is.na(x)) {return('NA')}
	else if(x == '0') {return('0')}
        else {return(omxQuotes(x))}
      })
      print(matrix, quote = FALSE)
   } else {
      print(matrix)
   }
} 

setMethod("t",  "MxSymmetricSparse",		
	function(x) {
		return(x)
	}
)

setMethod("nnzero", "MxSymmetricSparse",
	function(x, na.counted = NA) {
	    duplicates <- sum(x@rowVector == x@colVector)	    
		return(length(x@rowVector) * 2 - duplicates)
	}
)
		
setMethod("as.matrix",  "MxSymmetricSparse",
	function(x, ...) {
		result <- matrix(0, nrow(x), ncol(x))
		result[x@rowVector + (x@colVector - 1)*nrow(x)] <- x@dataVector
		result[x@colVector + (x@rowVector - 1)*ncol(x)] <- x@dataVector		
		return(result)
	}
)

setMethod("nrow", "MxSymmetricSparse",
	function(x) {
	    return(x@nrow)
	}
)

setMethod("ncol", "MxSymmetricSparse",
	function(x) {
		return(x@ncol)
	}
)



setMethod("[", "MxSymmetricSparse",
	function(x, i, j, ..., drop = FALSE) {
       	if (is.character(i)) {
       		if(length(x@dimNames) == 0) {
       			stop("no 'dimnames' attribute for array")
       		} else {
       			rownames <- x@dimNames[[1]]
       			i <- match(i, rownames)	
       		}
       	}
       	if (is.character(j)) {
       		if(length(x@dimNames) == 0) {       				stop("no 'dimnames' attribute for array")
       		} else {
       			colnames <- x@dimNames[[2]]
       			j <- match(j, colnames)	
       		}
       	}        							
	    if (i > x@nrow || j > x@ncol) {
	        stop("subscript out of bounds")   	
		} else if (length(x@rowVector) == 0) {
			return(0)
        } else {
        	if (i > j) {
        		tmp <- i
			    i <- j
			    j <- tmp        	        	
        	}
			ubound <- findInterval(i, x@rowVector)
			if ((ubound > 0) && x@rowVector[ubound] == i) {
				lbound <- findInterval(i - 1, x@rowVector) + 1
				index <- findInterval(j, x@colVector[lbound:ubound])
				offset <- index + lbound - 1
				if (offset > 0 && x@colVector[offset] == j) {
					return(x@dataVector[offset])
				}
			} 
		}
		return(0)
    }
)

setReplaceMethod("[", "MxSymmetricSparse", 
	function(x, i, j, value) {
       	if (is.character(i)) {
       		if(length(x@dimNames) == 0) {
       			stop("no 'dimnames' attribute for array")
       		} else {
       			rownames <- x@dimNames[[1]]
       			i <- match(i, rownames)	
       		}
       	}
       	if (is.character(j)) {
       		if(length(x@dimNames) == 0) {       				stop("no 'dimnames' attribute for array")
       		} else {
       			colnames <- x@dimNames[[2]]
       			j <- match(j, colnames)	
       		}
       	}		
	    if (i > x@nrow || j > x@ncol) {
	        stop("subscript out of bounds")   	
	    } else if (i > j) {
        		tmp <- i
			    i <- j
			    j <- tmp        	        	
        }	
		if (length(x@rowVector) == 0) {
			x@rowVector <- c(i)
			x@colVector <- c(j)
			x@dataVector <- c(value)
        } else {
			ubound <- findInterval(i, x@rowVector)
			if ((ubound > 0) && x@rowVector[ubound] == i) {
				lbound <- findInterval(i - 1, x@rowVector) + 1
				index <- findInterval(j, x@colVector[lbound:ubound])
				offset <- index + lbound - 1
				if (offset > 0 && x@colVector[offset] == j) {
				    if (!is.na(value) && value == 0) {
					   x@rowVector <- x@rowVector[-offset]
				       x@colVector <- x@colVector[-offset]
				       x@dataVector <- x@dataVector[-offset]
				    } else {
  					   x@rowVector[offset] <- i
					   x@colVector[offset] <- j
					   x@dataVector[offset] <- value
				    }
				} else if (is.na(value) || value != 0) {
					x@rowVector <- append(x@rowVector, i, after = offset)
					x@colVector <- append(x@colVector, j, after = offset)
					x@dataVector <- append(x@dataVector, value, after = offset)
				}
			} else if (is.na(value) || value != 0) {
				x@rowVector <- append(x@rowVector, i, after = ubound)
				x@colVector <- append(x@colVector, j, after = ubound)
				x@dataVector <- append(x@dataVector, value, after = ubound)
			}
		}
		return(x)
    }
)
