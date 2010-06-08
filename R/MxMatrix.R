#
#   Copyright 2007-2010 The OpenMx Project
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


nnzero <- function(matrix) {
	return(length(which(matrix != 0)))
}

single.na <- function(a) {
	return((length(a) == 1) && 
		(is.list(a) || is.vector(a) || is.matrix(a)) && 
		(is.na(a) == TRUE))
}

all.na <- function(a) {
	return((length(a) > 0) &&
		(is.list(a) || is.vector(a) || is.matrix(a)) && 
		(all(sapply(a, is.na))))	
}

verifySquare <- function(.Object) {
	if (nrow(.Object@labels) != ncol(.Object@labels)) { 
		stop(paste("Labels matrix of MxMatrix", 
				omxQuotes(.Object@name), "is not square"), call.=FALSE)
	}
	if (nrow(.Object@values) != ncol(.Object@values)) {
		stop(paste("Values matrix of MxMatrix", 
				omxQuotes(.Object@name), "is not square"), call.=FALSE)
	}
	if (nrow(.Object@free) != ncol(.Object@free)) { 
		stop(paste("Free matrix of MxMatrix", 
				omxQuotes(.Object@name), "is not square"), call.=FALSE)
	}
}

setGeneric("omxVerifyMatrix", function(.Object) { 
	return(standardGeneric("omxVerifyMatrix")) 
} )

setGeneric("omxSymmetricMatrix", function(.Object) {
	return(standardGeneric("omxSymmetricMatrix"))
})

setGeneric("omxSquareMatrix", function(.Object) {
	return(standardGeneric("omxSquareMatrix"))
})

#
# The Matrix package returns a non-symmetric matrix
# when you modify a symmetric matrix.
# This prevents us from using symmetric matrices from the package.
#
setClass(Class = "MxMatrix",
	representation = representation(
		labels = "matrix", values = "matrix", 
		free = "matrix", name = "character", 
		lbound = "matrix", ubound = "matrix", "VIRTUAL"))		
		
setMethod("initialize", "MxMatrix",
	function(.Object, labels, values, free, lbound, ubound, name) {
		.Object@labels = labels
		.Object@values = values
		.Object@free = free
		.Object@name = name
		.Object@lbound = lbound
		.Object@ubound = ubound
		return(.Object)
	}
)		

setMethod("omxVerifyMatrix", "MxMatrix",
	function(.Object) {
		if (!all(dim(.Object@labels) == dim(.Object@values))) {
			stop(paste("Labels and values matrices of", 
				omxQuotes(.Object@name), 
				"have different dimensions"), call.=FALSE)
		}
		if (!all(dim(.Object@labels) == dim(.Object@free))) {
			stop(paste("Labels and free matrices of", 
				omxQuotes(.Object@name), 
				"have different dimensions"), call.=FALSE)
		}
		if (omxSquareMatrix(.Object)) {
			verifySquare(.Object)
		}
	}
)

setMethod("omxSymmetricMatrix", "MxMatrix",
	function(.Object) { return(FALSE) }
)

setMethod("omxSquareMatrix", "MxMatrix",
	function(.Object) { return(FALSE) }
)

setMethod("nrow", "MxMatrix",
	function(x) {
	    return(nrow(x@values))
	}
)

setMethod("ncol", "MxMatrix",
	function(x) {
	    return(ncol(x@values))
	}
)

setMethod("length", "MxMatrix",
	function(x) {
	    return(nrow(x) * ncol(x))
	}
)

setMethod("[", "MxMatrix",
	function(x, i, j, ..., drop = FALSE) {
		if(!is.null(match.call()$name)) {
			name <- match.call()$name
		} else {
			name <- x@name
		}
		labels <- as.matrix(x@labels[i,j])
		values <- as.matrix(x@values[i,j])		
		free <- as.matrix(x@free[i,j])
		lbound <- as.matrix(x@lbound[i,j])
		ubound <- as.matrix(x@ubound[i,j])
		if(!missing(i) && missing(j) && length(i) == 1) {
			labels <- t(labels)
			values <- t(values)
			free <- t(free)
			lbound <- t(lbound)
			ubound <- t(ubound)
		}
		type <- class(x)[[1]]
		nrow <- nrow(labels)
		ncol <- ncol(labels)
		newMatrix <- tryCatch(suppressWarnings(
			new(type, name, values, free, labels, 
				lbound, ubound, nrow, ncol, FALSE)),
				error = function(e) new("FullMatrix", name, 
					values, free, labels, lbound, ubound, 
					nrow, ncol, FALSE))
		return(newMatrix)
	}
)

setReplaceMethod("[", "MxMatrix", 
	function(x, i, j, value) {
		if(!is(value,"MxMatrix")) {
			stop("right-hand side must be MxMatrix object")
		}
		x@values[i,j] <- value@values
		x@free[i,j] <- value@free
		x@labels[i,j] <- value@labels
		x@lbound[i,j] <- value@lbound
		x@ubound[i,j] <- value@ubound
		return(x)
	}
)

setMethod("dimnames", "MxMatrix",
	function(x) { return(dimnames(x@values)) }
)


setReplaceMethod("dimnames", "MxMatrix",
	function(x, value) {
		dimnames(x@values) <- value
		dimnames(x@free) <- value
		dimnames(x@labels) <- value
		dimnames(x@lbound) <- value
		dimnames(x@ubound) <- value
		return(x)
	}
)

matrixTypes <- c("Diag", "Full", "Iden", "Lower", "Stand", "Sdiag", "Symm", "Unit", "Zero")
squareMatrices <- c("Diag", "Iden", "Lower", "Stand", "Sdiag", "Symm")


matrixCheckDims <- function(type, values, free, labels, lbound, ubound, nrow, ncol) {
	if (is.matrix(values)) {
		nrow <- nrow(values)
		ncol <- ncol(values)
	} else if (is.matrix(labels)) {
		nrow <- nrow(labels)
		ncol <- ncol(labels)
	} else if (is.matrix(free)) {
		nrow <- nrow(free)
		ncol <- ncol(free)	
	} else if (is.matrix(lbound)) {
		nrow <- nrow(lbound)
		ncol <- ncol(lbound)
	} else if (is.matrix(ubound)) {
		nrow <- nrow(ubound)
		ncol <- ncol(ubound)
	} else if ((type == "Full") && ((is.na(nrow) && !is.na(ncol)) ||
		(!is.na(nrow) && is.na(ncol)))) {
		allLengths <- c(length(values), length(labels), 
			length(free), length(lbound), length(ubound))
		nonSingle <- allLengths[allLengths > 1]
		if (length(nonSingle) > 0) {
			first <- nonSingle[[1]]
			if (all(sapply(nonSingle, function(x) { x == first }))) {
				if (is.na(nrow)) {
					nrow <- ceiling(first / ncol)
				} else {
					ncol <- ceiling(first / nrow)				
				}
			}		
		}
	}
	return(c(nrow, ncol))
}

mxMatrix <- function(type = "Full", nrow = NA, ncol = NA, 
	free = FALSE, values = NA, labels = NA, 
	lbound = NA, ubound = NA, byrow = getOption('mxByrow'), 
	dimnames = NA, name = NA) {
	if (all.na(values)) { values <- as.numeric(values) }
	if (all.na(labels)) { labels <- as.character(labels) }
	if (all.na(lbound)) { lbound <- as.numeric(lbound) }
	if (all.na(ubound)) { ubound <- as.numeric(ubound) }
	if (single.na(name)) {
		name <- omxUntitledName()
	}
	omxVerifyName(name, 0)
	if (!is.character(name)) {
		stop(paste("'name' argument must",
			"be a character vector in", 
			deparse(width.cutoff = 400L, sys.call())), call. = FALSE)
	}
	matrixCheckErrors(type, values, free, labels, lbound, ubound, nrow, ncol, name)
	checkDims <- matrixCheckDims(type, values, free, labels, lbound, ubound, nrow, ncol)
	nrow <- checkDims[[1]]
	ncol <- checkDims[[2]]
	if (type %in% squareMatrices) {
		if (is.na(nrow) && is.na(ncol)) {
			stop(paste("Either 'nrow' or 'ncol'",
				"must be specified on a",
				"square matrix in", 
				deparse(width.cutoff = 400L, sys.call())), call. = FALSE)
		} else if (is.na(nrow)) {
			nrow <- ncol
		} else if (is.na(ncol)) {
			ncol <- nrow
		}
	} else if (is.na(nrow) || is.na(ncol)) {
		stop(paste("Both 'nrow' and 'ncol'",
					"must be specified on a",
					"non-square matrix in",
					deparse(width.cutoff = 400L, sys.call())), call. = FALSE)
	}
	values <- as.numeric.preserve(values)
	lbound <- as.numeric.preserve(lbound)
	ubound <- as.numeric.preserve(ubound)
	typeName <- paste(type, "Matrix", sep="")
	newMatrix <- new(typeName, name, values, free, labels, 
			lbound, ubound, nrow, ncol, byrow)
	if(length(dimnames) == 1 && is.na(dimnames)) {
	} else {
		dimnames(newMatrix) <- dimnames
	}
	omxVerifyMatrix(newMatrix)
	return(newMatrix)
}

as.numeric.preserve <- function(x, ...) {
	if (is.matrix(x)) {
		return(matrix(as.numeric(x, ...), nrow(x), ncol(x)))
	} else {
		return(as.numeric(x, ...))
	}
}

matrixCheckArgument <- function(arg, name) {
	if (is.list(arg) || isS4(arg)) {
		stop(paste(omxQuotes(name), "argument to mxMatrix",
			"must be a scalar, a vector, or a matrix in",
			deparse(width.cutoff = 400L, sys.call(-2))), call. = FALSE)
	}
}


matrixCheckErrors <- function(type, values, free, labels, lbound, ubound, nrow, ncol, name) {
	if (is.na(match(type, matrixTypes))) {
		stop(paste("Type must be one of:", 
			paste(matrixTypes, collapse=" "),
			"in", deparse(width.cutoff = 400L, sys.call(-1))), call. = FALSE)
	}
	matrixCheckArgument(values, 'values')
	matrixCheckArgument(free, 'free')
	matrixCheckArgument(labels, 'labels')
	matrixCheckArgument(lbound, 'lbound')
	matrixCheckArgument(ubound, 'ubound')
	if (!is.numeric(values)) {
		stop(paste("'values' argument to mxMatrix",
			"must be of numeric type in", 
			deparse(width.cutoff = 400L, sys.call(-1))), call. = FALSE)
	}
	if (!is.logical(free)) {
		stop(paste("'free' argument to mxMatrix",
			"must be of logical type in", 
			deparse(width.cutoff = 400L, sys.call(-1))), call. = FALSE)
	}
	if (!is.character(labels)) {
		stop(paste("'labels' argument to mxMatrix",
			"must be of character type in", 
			deparse(width.cutoff = 400L, sys.call(-1))), call. = FALSE)
	}
	if (!is.numeric(lbound)) {
		stop(paste("'lbound' argument to mxMatrix",
			"must be of numeric type in", 
			deparse(width.cutoff = 400L, sys.call(-1))), call. = FALSE)
	}
	if (!is.numeric(ubound)) {
		stop(paste("'ubound' argument to mxMatrix",
			"must be of numeric type in", 
			deparse(width.cutoff = 400L, sys.call(-1))), call. = FALSE)
	}
	if (!(is.na(nrow) || (is.numeric(nrow) && length(nrow) == 1))) {
		stop(paste("'nrow' argument to mxMatrix",
			"must be either NA or a single value in", 
			deparse(width.cutoff = 400L, sys.call(-1))), call. = FALSE)
	}
	if (!(is.na(ncol) || (is.numeric(ncol) && length(ncol) == 1))) {
		stop(paste("'ncol' argument to mxMatrix",
			"must be either NA or a single value in",
			deparse(width.cutoff = 400L, sys.call(-1))), call. = FALSE)
	}
	if ((is.matrix(values) || is.matrix(free) || is.matrix(labels) 
		|| is.matrix(lbound) || is.matrix(ubound)) &&
		(!is.na(nrow))) {
		warning(paste("'nrow' is disregarded",
			"for mxMatrix constructor in", 
			deparse(width.cutoff = 400L, sys.call(-1))), call. = FALSE)
	}
	if ((is.matrix(values) || is.matrix(free) || is.matrix(labels)
		|| is.matrix(lbound) || is.matrix(ubound)) &&
		(!is.na(ncol))) {
		warning(paste("'ncol' is disregarded",
			"for mxMatrix constructor in", 
			deparse(width.cutoff = 400L, sys.call(-1))), call. = FALSE)
	}
	dimensions <- sapply(list(values, free, labels, lbound, ubound), dim)
	dimensions <- dimensions[!sapply(dimensions, is.null)]
	if (length(dimensions) > 1) {
		allEqual <- sapply(dimensions, function(x) { x == dimensions[[1]] })
		if(!all(allEqual)) {
			stop(paste("Two matrices provided",
				"to mxMatrix are not of identical",
				"dimensions in", deparse(width.cutoff = 400L, sys.call(-1))), call. = FALSE)
		}
	}
	lapply(labels, omxVerifyReference, -2)
	if(any(is.na(free))) {
		stop(paste("'free' argument to mxMatrix",
			"cannot contain an NA in",
			deparse(width.cutoff = 400L, sys.call(-1))), call. = FALSE)
	}
}

matrixParameters <- function(free, labels, lbound, ubound, 
	result, matrixNumber, isSymmetric) {
	if (all(free == FALSE)) {
		return(result)
	}
	if (isSymmetric) {
		triangle <- upper.tri(free, diag=TRUE)
		select <- free & triangle
	} else {
		select <- free
	}
	parameterNames <- labels[select]
	rows <- row(labels)[select]
	cols <- col(labels)[select]
	lbound <- lbound[select]
	ubound <- ubound[select]
	for(i in 1:length(parameterNames)) {
		parameterName <- parameterNames[i]
		row <- rows[i] - 1L
		col <- cols[i] - 1L
		minBounds <- lbound[i]
		maxBounds <- ubound[i]
		if (is.na(parameterName)) {
			result[[length(result)+1]] <-  list(minBounds, maxBounds, 
				c(matrixNumber, row, col))
			if (isSymmetric && row != col) {
				original <- result[[length(result)]]
				original[[length(original) + 1]] <- c(matrixNumber, col, row)
				result[[length(result)]] <- original
			}
			names(result)[[length(result)]] <- as.character(NA)
		} else if (length(grep(omxSeparatorChar, parameterName, fixed = TRUE)) == 0) {
			if (!is.null(result[[parameterName]])) {
				original <- result[[parameterName]]
				original[[length(original) + 1]] <- c(matrixNumber, row, col)
				if (isSymmetric && row != col) {
					original[[length(original) + 1]] <- c(matrixNumber, col, row)
				}
				result[[parameterName]] <- original
			} else {
				result[[parameterName]] <- list('min' = minBounds, 'max' = maxBounds, 
					c(matrixNumber, row, col))
				if (isSymmetric && row != col) {
					original <- result[[parameterName]]
					original[[length(original) + 1]] <- c(matrixNumber, col, row)
					result[[parameterName]] <- original
				}
			}
		}
	}
	return(result)
}

# Definition variables is a list:
# each entry of the list is a sublist of length 3 or more
# first entry of the sublist: data number
# second entry of the sublist: column number
# remaining entries of the sublist: c(matrix number, row, column)
matrixDefinitions <- function(free, labels, result, defLocations, matrixNumber) {
	select <- !is.na(labels)
	if (all(select == FALSE)) {
		return(result)
	}
	defNames <- names(defLocations)
	parameterNames <- labels[select]
	rows <- row(labels)[select]
	cols <- col(labels)[select]
	for(i in 1:length(parameterNames)) {
		parameterName <- parameterNames[i]
		row <- rows[i] - 1L
		col <- cols[i] - 1L
		if (parameterName %in% defNames) {
			if (!is.null(result[[parameterName]])) {
				original <- result[[parameterName]]
				original[[length(original) + 1]] <- c(matrixNumber, row, col)
				result[[parameterName]] <- original
			} else {
				dataNumber <- as.integer(defLocations[[parameterName]][[1]])
				columnNumber <- as.integer(defLocations[[parameterName]][[2]])
				result[[parameterName]] <- list(dataNumber, columnNumber, 
						c(matrixNumber, row, col))
			}
		}
	}
	return(result)
}

generateMatrixValuesHelper <- function(mxMatrix) {
	return(mxMatrix@values)
}

generateParameterListHelper <- function(mxMatrix,
	result, matrixNumber) {
	free <- mxMatrix@free
	labels <- mxMatrix@labels
	lbound <- mxMatrix@lbound
	ubound <- mxMatrix@ubound
	isSymmetric <- omxSymmetricMatrix(mxMatrix)
	result <- matrixParameters(free, labels, lbound,
		ubound, result, matrixNumber, isSymmetric)
	return(result)
}

generateDefinitionListHelper <- function(mxMatrix,
	result, defLocations, matrixNumber) {
	free <- mxMatrix@free
	labels <- mxMatrix@labels
	result <- matrixDefinitions(free, labels,
		result, defLocations, matrixNumber)
	return(result)
}

displayMatrix <- function(mxMatrix) {
	type <- class(mxMatrix)[[1]]
	cat(type, omxQuotes(mxMatrix@name), '\n')
	cat("\n")
	nolabels <- all(is.na(mxMatrix@labels))
	if(nolabels == FALSE) {
		cat("@labels\n")
		print(mxMatrix@labels)
		cat("\n")
	} else {
		cat("@labels: No labels assigned.\n\n")
	}
	cat("@values\n")
	print(mxMatrix@values)
	cat("\n")
	noFree <- all(mxMatrix@free == FALSE)
	if(noFree == FALSE) {
		cat("@free\n")
		print(mxMatrix@free)
		cat("\n")
	} else {
		cat("@free: No free parameters.\n\n")
	}
	if(!all(is.na(mxMatrix@lbound))) {
		cat("@lbound\n")
		print(mxMatrix@lbound)
		cat("\n")
	} else {
		cat("@lbound: No lower bounds assigned.\n\n")
	}
	if(!all(is.na(mxMatrix@ubound))) {
		cat("@ubound\n")
		print(mxMatrix@ubound)
		cat("\n")
	} else {
		cat("@ubound: No upper bounds assigned.\n\n")
	}
}

setMethod("print", "MxMatrix", function(x,...) { displayMatrix(x) })
setMethod("show", "MxMatrix", function(object) { displayMatrix(object) })
