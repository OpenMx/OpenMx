single.na <- function(a) {
	return((length(a) == 1) && 
		(is.list(a) || is.vector(a)) && 
		(is.na(a) == TRUE))
}

is.Matrix <- function(a) {
	return(is.matrix(a) || is(a, "Matrix"))	
}

verifySquare <- function(.Object) {
	if (nrow(.Object@specification) != ncol(.Object@specification)) { 
		stop(paste("Specification matrix of MxMatrix", 
				omxQuotes(.Object@name), "is not square"), call.=FALSE)
	}
	if (nrow(.Object@values) != ncol(.Object@values)) {
		stop(paste("Values matrix of MxMatrix", 
				omxQuotes(.Object@name), "is not square"), call.=FALSE)
	}		
}

setGeneric("omxVerifyMatrix", function(.Object) { 
	return(standardGeneric("omxVerifyMatrix")) 
} )

#
# The Matrix package returns a non-symmetric matrix
# when you modify a symmetric matrix.
# This prevents us from using symmetric matrices from the package.
#
setClass(Class = "MxSymmetricMatrix",
	representation = representation(
		specification = "MxSymmetricSparse",
		values = "MxSymmetricSparse", name = "character", "VIRTUAL"))

setClass(Class = "MxNonSymmetricMatrix",
	representation = representation(
		specification = "MxSparseMatrix",
		values = "Matrix", name = "character", "VIRTUAL"))
		
setClassUnion("MxMatrix",
    c("MxSymmetricMatrix", "MxNonSymmetricMatrix"))
		
		
setMethod("initialize", "MxMatrix",
	function(.Object, specification, values, name) {
		.Object@specification = specification
		.Object@values = values
		.Object@name = name
		return(.Object)
	}
)		

setMethod("omxVerifyMatrix", "MxMatrix",
	function(.Object) {
		if (nrow(.Object@specification) != nrow(.Object@values)) {
			stop(paste("Specification and values matrices of", 
				omxQuotes(.Object@name), "have different dimensions"), call.=FALSE)
		}
		if (ncol(.Object@specification) != ncol(.Object@values)) {
			stop(paste("Specification and values matrices of", 
				omxQuotes(.Object@name), "have different dimensions"), call.=FALSE)
		}		
	}
)

setMethod("nrow", "MxMatrix",
	function(x) {
	    return(nrow(x@specification))
	}
)

setMethod("ncol", "MxMatrix",
	function(x) {
	    return(ncol(x@specification))
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

matrixTypes <- c("Diag", "Full", "Iden", "Symm", "Unit", "Zero")
squareMatrices <- c("Diag", "Iden", "Symm")


mxMatrix <- function(type = "Full", values = NA, 
	specification = NA, nrow = NA, ncol = NA,
	byrow = FALSE, free = FALSE, name = NA) {
	if (is.na(match(type, matrixTypes))) {
		stop(paste("Type must be one of:", paste(matrixTypes, collapse=" ")))
	}
	if ((is.Matrix(values) || is.Matrix(specification)) &&
		(!is.null(match.call()$nrow))) {
		warning("\'nrow\' is disregarded for mxMatrix constructor")
	}
	if ((is.Matrix(values) || is.Matrix(specification)) &&
		(!is.null(match.call()$ncol))) {
		warning("\'ncol\' is disregarded for mxMatrix constructor")
	}
	if (is.Matrix(values) && is.Matrix(specification) &&
		!all(dim(values) == dim(specification))) {
		stop("Values and specification matrices are not of identical dimensions")
	}
	if (is.Matrix(values)) {
		nrow <- nrow(values)
		ncol <- ncol(values)
	} else if (is.Matrix(specification)) {
		nrow <- nrow(specification)
		ncol <- ncol(specification)
	}
	if (type %in% squareMatrices) {
		if (is.na(nrow) && is.na(ncol)) {
			stop("Either nrow or ncol must be specified on a square matrix")
		} else if (is.na(nrow)) {
			nrow <- ncol
		} else if (is.na(ncol)) {
			ncol <- nrow
		}
	} else if (is.na(nrow) || is.na(ncol)) {
		stop("Both nrow and ncol must be specified on a non-square matrix")
	}
	if (byrow) {
		values <- omxCreateByRow(values, nrow, ncol, type)
		specification <- omxCreateByRow(specification, nrow, ncol, type)
	}
	if (is.na(name)) {
		name <- omxUntitledName()
	}
	if (!is.character(name)) {
		stop("\'name\' must be a character vector!")
	}
	typeName <- paste(type, "Matrix", sep="")
	return(new(typeName, name, values, specification, 
			nrow, ncol, byrow, free))
}

omxCreateByRow <- function(source, nrow, ncol, type) {
	if (!is.na(source) && is.vector(source)) {
		if (type %in% squareMatrices) {
			res <- matrix(0, nrow, ncol)
			res[upper.tri(res, diag=TRUE)] <- source
			res[lower.tri(res)] <- res[upper.tri(res)]
			return(res)
		} else {
			return(matrix(source, nrow, ncol, TRUE))
		}
	} else {
		return(source)
	}
}

omxSparseMatrixParameters <- function(specification, bounds, result, 
	defNames, matrixNumber, reverse=FALSE) {
	if (length(specification@dataVector) == 0) {
		return(result)
	}
	for(i in 1:length(specification@dataVector)) {
	    if (reverse == FALSE || specification@rowVector[i] != 
	    	specification@colVector[i]) {
		    parameterName <- as.character(specification@dataVector[[i]])
		    if (reverse) {
			    col <- specification@rowVector[i] - 1
			    row <- specification@colVector[i] - 1
			} else {
			    row <- specification@rowVector[i] - 1
			    col <- specification@colVector[i] - 1			
			}
			boundsLookup <- omxLocateBounds(bounds, parameterName)
			minBounds <- boundsLookup[[1]]
			maxBounds <- boundsLookup[[2]]
			if (is.na(parameterName)) {
				result[[length(result)+1]] <-  list(minBounds, maxBounds, 
					c(matrixNumber, row, col))
			} else if (!(parameterName %in% defNames)) {
				if (!is.null(result[[parameterName]])) {
					original <- result[[parameterName]]
					original[[length(original) + 1]] <- c(matrixNumber, row, col)
					result[[parameterName]] <- original
				} else {
					result[[parameterName]] <- list(minBounds, maxBounds, 
						c(matrixNumber, row,col))
				}
			}
		}
	}
	return(result)
}

omxSparseMatrixDefinitions <- function(specification, bounds, result, 
	defLocations, matrixNumber, reverse=FALSE) {
	if (length(specification@dataVector) == 0) {
		return(result)
	}
	defNames <- names(defLocations)
	for(i in 1:length(specification@dataVector)) {
	    if (reverse == FALSE || specification@rowVector[i] != 
	    	specification@colVector[i]) {
		    parameterName <- as.character(specification@dataVector[[i]])
		    if (reverse) {
			    col <- specification@rowVector[i] - 1
			    row <- specification@colVector[i] - 1
			} else {
			    row <- specification@rowVector[i] - 1
			    col <- specification@colVector[i] - 1			
			}
			if (parameterName %in% defNames) {
				if (!is.null(result[[parameterName]])) {
					original <- result[[parameterName]]
					dataNumber <- 
					original[[length(original) + 1]] <- c(matrixNumber, row, col)
					result[[parameterName]] <- original
				} else {
					dataNumber <- defLocations[[parameterName]][[1]]
					columnNumber <- defLocations[[parameterName]][[2]]
					result[[parameterName]] <- list(dataNumber, columnNumber, 
						c(matrixNumber, row,col))
				}
			}
		}
	}
	return(result)
}

generateMatrixListHelper <- function(mxMatrix) {
	values <- mxMatrix@values
    if(is(values,"MxSymmetricSparse")) {
		return(as.matrix(values))    #return(Matrix(as.matrix(values)))    
    } else {
		return(as.matrix(values))
    }
}

omxGenerateParameterListHelper <- function(mxMatrix, bounds,
	result, defNames, matrixNumber) {
	specification <- mxMatrix@specification
	result <- omxSparseMatrixParameters(specification, bounds,
		result, defNames, matrixNumber)
	if(is(specification,"MxSymmetricSparse")) {
		result <- omxSparseMatrixParameters(specification, bounds, 
			result, defNames, matrixNumber, TRUE)
	}
	return(result)
}

omxGenerateDefinitionListHelper <- function(mxMatrix, bounds,
	result, defLocations, matrixNumber) {
	specification <- mxMatrix@specification
	result <- omxSparseMatrixDefinitions(specification, bounds,
		result, defLocations, matrixNumber)
	if(is(specification,"MxSymmetricSparse")) {
		result <- omxSparseMatrixDefinitions(specification, bounds, 
			result, defLocations, matrixNumber, TRUE)
	}
	return(result)
}

omxDisplayMatrix <- function(mxMatrix) {
   cat("MxMatrix", omxQuotes(mxMatrix@name), '\n')
   cat("\n")
   cat("Specification matrix:\n")
   print(mxMatrix@specification, use.quotes = TRUE)
   cat("\n")
   cat("Values matrix:\n")
   values <- mxMatrix@values
   if(is(values, "sparseMatrix")) {
      print(as(values, 'matrix'))
   } else if(is(values, "MxSymmetricSparse")) {
      print(values, use.quotes = FALSE)
   } else {
      print(values)
   }
   cat("\n")
}

setMethod("print", "MxMatrix", function(x,...) { omxDisplayMatrix(x) })
setMethod("show", "MxMatrix", function(object) { omxDisplayMatrix(object) })
