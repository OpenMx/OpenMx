#
#   Copyright 2007-2009 The OpenMx Project
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
		(is.list(a) || is.vector(a)) && 
		(is.na(a) == TRUE))
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

#
# The Matrix package returns a non-symmetric matrix
# when you modify a symmetric matrix.
# This prevents us from using symmetric matrices from the package.
#
setClass(Class = "MxSymmetricMatrix",
	representation = representation(
		labels = "matrix", values = "matrix", 
		free = "matrix", name = "character", 
		lbound = "matrix", ubound = "matrix", "VIRTUAL"))

setClass(Class = "MxNonSymmetricMatrix",
	representation = representation(
		labels = "matrix", values = "matrix", 
		free = "matrix", name = "character", 
		lbound = "matrix", ubound = "matrix", "VIRTUAL"))
		
setClassUnion("MxMatrix",
    c("MxSymmetricMatrix", "MxNonSymmetricMatrix"))
		
		
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
		if(!all(dim(.Object@labels) == dim(.Object@values))) {
			stop(paste("Labels and values matrices of", 
				omxQuotes(.Object@name), 
				"have different dimensions"), call.=FALSE)
		}
		if(!all(dim(.Object@labels) == dim(.Object@free))) {
			stop(paste("Labels and free matrices of", 
				omxQuotes(.Object@name), 
				"have different dimensions"), call.=FALSE)
		}
	}
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

setMethod("[", "MxMatrix",
	function(x, i, j, ..., drop = FALSE) {
		return(x@values[i,j])
	}
)

setReplaceMethod("[", "MxMatrix", 
	function(x, i, j, value) {
		if(is(value, "numeric")) {
			x@values[i,j] <- value
		} else if(is(value, "logical")) {
			x@free[i,j] <- value
		} else if(is(value, "character")) {
			x@labels[i,j] <- value
		} else {
			stop(paste("Unknown type", type(value), "for MxMatrix",
				omxQuotes(x@name), "element replacement."), call. = FALSE)
		}
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

matrixTypes <- c("Diag", "Full", "Iden", "Symm", "Unit", "Zero")
squareMatrices <- c("Diag", "Iden", "Symm")


mxMatrix <- function(type = "Full", values = NA, free = FALSE, 
	labels = NA, lbound = NA, ubound = NA, nrow = NA, ncol = NA, 
	byrow = FALSE, dimnames = NA, name = NA) {
	matrixCheckErrors(type, values, free, labels, lbound, ubound, nrow, ncol)
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
	if (is.na(name)) {
		name <- omxUntitledName()
	}
	omxVerifyName(name)
	if (!is.character(name)) {
		stop("\'name\' must be a character vector!")
	}
	fiveMatrices <- convertVFN(values, free, labels, lbound, ubound, nrow, ncol)
	values <- fiveMatrices[[1]]
	free <- fiveMatrices[[2]]
	labels <- fiveMatrices[[3]]
	lbound <- fiveMatrices[[4]]
	ubound <- fiveMatrices[[5]]
	typeName <- paste(type, "Matrix", sep="")
	newMatrix <- new(typeName, name, values, free, labels, 
			lbound, ubound, nrow, ncol, byrow)
	if(length(dimnames) == 1 && is.na(dimnames)) {
	} else {
		dimnames(newMatrix) <- dimnames
	}
	return(newMatrix)
}

matrixCheckErrors <- function(type, values, free, labels, lbound, ubound, nrow, ncol) {
	if (is.na(match(type, matrixTypes))) {
		stop(paste("Type must be one of:", 
			paste(matrixTypes, collapse=" ")), call. = FALSE)
	}
	if ((is.matrix(values) || is.matrix(free) || is.matrix(labels) 
		|| is.matrix(lbound) || is.matrix(ubound)) &&
		(!is.na(nrow))) {
		warning("\'nrow\' is disregarded for mxMatrix constructor")
	}
	if ((is.matrix(values) || is.matrix(free) || is.matrix(labels)
		|| is.matrix(lbound) || is.matrix(ubound)) &&
		(!is.na(ncol))) {
		warning("\'ncol\' is disregarded for mxMatrix constructor")
	}
	dimensions <- sapply(list(values, free, labels, lbound, ubound), dim)
	dimensions <- dimensions[!sapply(dimensions, is.null)]
	if (length(dimensions) > 1) {
		allEqual <- sapply(dimensions, function(x) { x == dimensions[[1]] })
		if(!all(allEqual)) {
			stop("Two matrices provided to the call are not of identical dimensions.", call. = FALSE)
		}
	}
}

convertVFN <- function(values, free, labels, lbound, ubound, nrow, ncol) {
	if (is.matrix(values)) {
		values <- matrix(as.numeric(values), nrow, ncol)
	} else if (is.vector(values)) {
		values <- as.numeric(values)
	} else {
		stop("\'values\' must be either a vector or a matrix", call. = FALSE)
	}
	if (is.matrix(free)) {
		free <- matrix(as.logical(free), nrow, ncol)
	} else if (is.vector(free)) {
		free <- as.logical(free)
	} else {
		stop("\'free\' must be either a vector or a matrix", call. = FALSE)
	}
	if (is.matrix(labels)) {
		labels <- matrix(as.character(labels), nrow, ncol)
	} else if (single.na(labels)) {
		labels <- as.character(NA)
	} else if (is.vector(labels)) {
		labels <- as.character(labels)
	} else {
		stop("\'labels\' must be either a vector or a matrix", call. = FALSE)
	}
	if (is.matrix(lbound)) {
		lbound <- matrix(as.numeric(lbound), nrow, ncol)
	} else if (is.vector(lbound)) {
		lbound <- as.numeric(lbound)
	} else {
		stop("\'lbound\' must be either a vector or a matrix", call. = FALSE)
	}
	if (is.matrix(ubound)) {
		ubound <- matrix(as.numeric(ubound), nrow, ncol)
	} else if (is.vector(ubound)) {
		ubound <- as.numeric(ubound)
	} else {
		stop("\'ubound\' must be either a vector or a matrix", call. = FALSE)
	}
	if(length(values) > 1 && any(is.na(values))) {
		stop("\'values\' cannot contain an NA", call. = FALSE)
	}
	if(any(is.na(free))) {
		stop("\'free\' cannot contain an NA", call. = FALSE)
	}
	return(list(values, free, labels, lbound, ubound))	
}


matrixParameters <- function(free, labels, lbound, ubound, 
	result, matrixNumber) {
	if (all(free == FALSE)) {
		return(result)
	}
	parameterNames <- labels[free]
	rows <- row(labels)[free]
	cols <- col(labels)[free]
	lbound <- lbound[free]
	ubound <- ubound[free]
	for(i in 1:length(parameterNames)) {
		parameterName <- parameterNames[i]
		row <- rows[i] - 1
		col <- cols[i] - 1
		minBounds <- lbound[i]
		maxBounds <- ubound[i]
		if (is.na(parameterName)) {
			result[[length(result)+1]] <-  list(minBounds, maxBounds, 
				c(matrixNumber, row, col))
			names(result)[[length(result)]] <- as.character(NA)
		} else if (length(grep(omxSeparatorChar, parameterName, fixed = TRUE)) == 0) {
			if (!is.null(result[[parameterName]])) {
				original <- result[[parameterName]]
				original[[length(original) + 1]] <- c(matrixNumber, row, col)
					result[[parameterName]] <- original
			} else {
				result[[parameterName]] <- list('min' = minBounds, 'max' = maxBounds, 
					c(matrixNumber, row,col))
			}
		}
	}
	return(result)
}

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
		row <- rows[i] - 1
		col <- cols[i] - 1
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
	result <- matrixParameters(free, labels, lbound,
		ubound, result, matrixNumber)
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
		cat("Labels matrix:\n")
		print(mxMatrix@labels)
		cat("\n")
	} else {
		cat("Labels matrix: No labels assigned.\n\n")
	}
	cat("Values matrix:\n")
	print(mxMatrix@values)
	cat("\n")
	noFree <- all(mxMatrix@free == FALSE)
	if(noFree == FALSE) {
		cat("Free matrix:\n")
		print(mxMatrix@free)
		cat("\n")
	} else {
		cat("Free matrix: No free parameters.\n\n")
	}
	if(!all(is.na(mxMatrix@lbound))) {
		cat("Lbound matrix:\n")
		print(mxMatrix@lbound)
		cat("\n")
	} else {
		cat("Lbound matrix: No lower bounds assigned.\n\n")
	}
	if(!all(is.na(mxMatrix@ubound))) {
		cat("Ubound matrix:\n")
		print(mxMatrix@ubound)
		cat("\n")
	} else {
		cat("Ubound matrix: No upper bounds assigned.\n\n")
	}
}

setMethod("print", "MxMatrix", function(x,...) { displayMatrix(x) })
setMethod("show", "MxMatrix", function(object) { displayMatrix(object) })
