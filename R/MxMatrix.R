verifySquare <- function(.Object) {
	if (nrow(.Object@specification) != ncol(.Object@specification)) 
		{ stop("Specification matrix is not square") }
	if (nrow(.Object@values) != ncol(.Object@values)) 
		{ stop("Values matrix is not square") }		
}

setGeneric("verify", function(.Object) { return(standardGeneric("verify")) } )

#
# The Matrix package returns a non-symmetric matrix
# when you modify a symmetric matrix.
# This prevents us from using symmetric matrices from the package.
#
setClass(Class = "MxSymmetricMatrix",
	representation = representation(
		specification = "MxSymmetricSparse",
		values = "MxSymmetricSparse", "VIRTUAL"))

setClass(Class = "MxNonSymmetricMatrix",
	representation = representation(
		specification = "MxSparseMatrix",
		values = "Matrix", "VIRTUAL"))
		
setClassUnion("MxMatrix",
    c("MxSymmetricMatrix", "MxNonSymmetricMatrix"))
		
		
setMethod("initialize", "MxMatrix",
	function(.Object, specification, values) {
		.Object@specification = specification
		.Object@values = values
		return(.Object)
	}
)		

setMethod("verify", "MxMatrix",
	function(.Object) {
		if (nrow(.Object@specification) != nrow(.Object@values)) 
			{ stop("Specification and values matrices have different dimensions") }
		if (ncol(.Object@specification) != ncol(.Object@values)) 
			{ stop("Specification and values matrices have different dimensions") }		
	}
)

setClass(Class = "ZeroMatrix",
	representation = representation(),
	contains = "MxNonSymmetricMatrix")

setMethod("initialize", "ZeroMatrix",
	function(.Object, nrow = 1, ncol = 1) {
		specification <- new("MxSparseMatrix", 0, nrow, ncol)
		values <- Matrix(0, nrow, ncol)
		return(callNextMethod(.Object, specification, values))
	}
)

setMethod("verify", "ZeroMatrix",
	function(.Object) {
		callNextMethod(.Object)	
		if(nnzero(.Object@specification) > 0) { stop("Specification matrix is not empty") } 
		if(nnzero(.Object@values) > 0) { stop("Values matrix is not empty") } 
	}
)

setClass(Class = "UnitMatrix",
	representation = representation(),
	contains = "MxNonSymmetricMatrix")
	
setMethod("initialize", "UnitMatrix",
	function(.Object, nrow = 1, ncol = 1) {
		specification <- new("MxSparseMatrix", 0, nrow, ncol)
		values <- Matrix(1, nrow, ncol)
		return(callNextMethod(.Object, specification, values))
	}
)

setMethod("verify", "UnitMatrix",
	function(.Object) {
		callNextMethod(.Object)		
		if(nnzero(.Object@specification) > 0) { stop("Specification matrix is not empty") } 
		if(nnzero(.Object@values - 1) > 0) { stop("Values matrix has non unit entries") } 
	}
)

setClass(Class = "IdenMatrix",
	representation = representation(),
	contains = "MxNonSymmetricMatrix")
	
setMethod("initialize", "IdenMatrix",
	function(.Object, nrow = 1) {
		specification <- new("MxSparseMatrix", 0, nrow, nrow)
		values <- Matrix(diag(nrow))
		return(callNextMethod(.Object, specification, values))
	}
)

setMethod("verify", "IdenMatrix",
	function(.Object) {
		callNextMethod(.Object)
		verifySquare(.Object)
		if(nnzero(.Object@specification) > 0) { stop("Specification matrix is not empty") } 		
		if(!suppressWarnings(all(.Object@values == diag(nrow(.Object@values))))) 
			{ stop("Values matrix is not the identity matrix") }
	}
)

setClass(Class = "DiagMatrix",
	representation = representation(),
	contains = "MxNonSymmetricMatrix")

setMethod("initialize", "DiagMatrix",
	function(.Object, data = 0, nrow = 1, free = FALSE) {
	    if (free) {
	    	specification <- diag(nrow)
	    	specification[specification == 1] <- NA
			specification <- new("MxSparseMatrix", specification)
		} else {
			specification <- new("MxSparseMatrix", 0, nrow, nrow)		
		}
		if (length(data) == 1) { values <- Matrix(data, nrow, nrow) }
   	 	else if (is.vector(data)) { values <- Matrix(data*diag(nrow)) }
    	else { values <- Matrix(data) }
		return(callNextMethod(.Object, specification, values))
	}
)

setMethod("verify", "DiagMatrix",
	function(.Object) {
		callNextMethod(.Object)
		verifySquare(.Object)
		values <- .Object@values
		if(suppressWarnings(nnzero(values - diag(diag(values,nrow=nrow(values),ncol=ncol(values))))) > 0)
			{ stop("Values matrix is not a diagonal matrix") }
	}
)

setClass(Class = "SymmMatrix",
	representation = representation(),
	contains = "MxSymmetricMatrix")
	
setMethod("initialize", "SymmMatrix",
	function(.Object, data = 0, nrow = 1, ncol = 1, free = FALSE) {
	    if (free) {
			specification <- new("MxSymmetricSparse", matrix(NA, nrow, ncol))
	    } else {
			specification <- new("MxSymmetricSparse", 0, nrow, ncol)
	    }
	    values <- new("MxSymmetricSparse", matrix(data,nrow = nrow, ncol = ncol))
		return(callNextMethod(.Object, specification, values))
	}
)


setClass(Class = "FullMatrix",
	representation = representation(),
	contains = "MxNonSymmetricMatrix")

setMethod("initialize", "FullMatrix",
	function(.Object, data = 0, nrow = 1, ncol = 1, free = FALSE) {
	    if (free) {
			specification <- new("MxSparseMatrix", matrix(NA, nrow, ncol))
	    } else {
			specification <- new("MxSparseMatrix", 0, nrow, ncol)
	    }
		values <- Matrix(data, nrow, ncol)
		return(callNextMethod(.Object, specification, values))
	}
)

setMethod("verify", "FullMatrix",
	function(.Object) {
		callNextMethod(.Object)
	}
)

setMethod("[", "MxMatrix",
	function(x, i, j, ..., drop = FALSE) {
		return(x@values[i,j])
    }
)

setReplaceMethod("[", "MxMatrix", 
	function(x, i, j, value) {
		x@values[i,j] <- value
		return(x)
    }
)


processSparseMatrix <- function(specification, result, matrixNumber, reverse=FALSE) {
	if (length(specification@dataVector) == 0) {
		return(result)
	}
	for(i in 1:length(specification@dataVector)) {
	    if (reverse == FALSE || specification@rowVector[i] != specification@colVector[i]) {
		    data <- as.character(specification@dataVector[[i]])
		    if (reverse) {
			    col <- specification@rowVector[i]
			    row <- specification@colVector[i]
			} else {
			    row <- specification@rowVector[i]
			    col <- specification@colVector[i]			
			}
			if (is.na(data)) {
				result[[length(result)+1]] <-  list(c(matrixNumber, row, col))
			} else {
				if (!is.null(result[[data]])) {
					original <- result[[data]]
					original[[length(original) + 1]] <- c(matrixNumber, row, col)
					result[[data]] <- original				
				} else {
					result[[data]] <- list(c(matrixNumber, row,col))
				}
			}
		}
	}
	return(result)
}

generateMatrixListHelper <- function(mxMatrix) {
	values <- mxMatrix@values
    if(is(values,"MxSymmetricSparse")) {
		return(Matrix(as.matrix(values)))    
    } else {
		return(values)    
    }
}

generaterParameterListHelper <- function(mxMatrix, result, matrixNumber) {
	specification <- mxMatrix@specification
	result <- processSparseMatrix(specification, result, matrixNumber)
	if(is(specification,"MxSymmetricSparse")) {
		result <- processSparseMatrix(specification, result, matrixNumber, TRUE)
	}
	return(result)
}
