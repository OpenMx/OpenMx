#
#   Copyright 2007-2014 The OpenMx Project
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

##' imxVerifyMatrix
##'
##' This is an internal function exported for those people who know
##' what they are doing.
##'
##' @param .Object .Object
##' @aliases
##' imxVerifyMatrix,DiagMatrix-method
##' imxVerifyMatrix,FullMatrix-method
##' imxVerifyMatrix,IdenMatrix-method
##' imxVerifyMatrix,LowerMatrix-method
##' imxVerifyMatrix,MxMatrix-method
##' imxVerifyMatrix,SdiagMatrix-method
##' imxVerifyMatrix,StandMatrix-method
##' imxVerifyMatrix,SymmMatrix-method
##' imxVerifyMatrix,UnitMatrix-method
##' imxVerifyMatrix,ZeroMatrix-method
setGeneric("imxVerifyMatrix", function(.Object) { 
	return(standardGeneric("imxVerifyMatrix")) 
} )

##' imxSymmetricMatrix
##' 
##' This is an internal function exported for those people who know
##' what they are doing.
##' 
##' @param .Object .Object
##' @aliases
##' imxSymmetricMatrix,LowerMatrix-method
##' imxSymmetricMatrix,MxMatrix-method
##' imxSymmetricMatrix,SdiagMatrix-method
##' imxSymmetricMatrix,StandMatrix-method
##' imxSymmetricMatrix,SymmMatrix-method
setGeneric("imxSymmetricMatrix", function(.Object) {
	return(standardGeneric("imxSymmetricMatrix"))
})

##' imxSquareMatrix
##'
##' This is an internal function exported for those people who know
##' what they are doing.
##'
##' @param .Object .Object
##' @aliases
##' imxSquareMatrix,DiagMatrix-method
##' imxSquareMatrix,IdenMatrix-method
##' imxSquareMatrix,LowerMatrix-method
##' imxSquareMatrix,MxMatrix-method
##' imxSquareMatrix,SdiagMatrix-method
##' imxSquareMatrix,StandMatrix-method
##' imxSquareMatrix,SymmMatrix-method
setGeneric("imxSquareMatrix", function(.Object) {
	return(standardGeneric("imxSquareMatrix"))
})

##' Create a matrix
##'
##' This is an internal function exported for those people who know
##' what they are doing.
##'
##' @param .Object the matrix
##' @param labels labels
##' @param values values
##' @param free free
##' @param lbound lbound
##' @param ubound ubound
##' @param nrow nrow
##' @param ncol ncol
##' @param byrow byrow
##' @param name name
##' @param ... Not used.
##' @aliases
##' imxCreateMatrix,MxMatrix-method
##' imxCreateMatrix,DiagMatrix-method
##' imxCreateMatrix,FullMatrix-method
##' imxCreateMatrix,IdenMatrix-method
##' imxCreateMatrix,LowerMatrix-method
##' imxCreateMatrix,SdiagMatrix-method
##' imxCreateMatrix,StandMatrix-method
##' imxCreateMatrix,SymmMatrix-method
##' imxCreateMatrix,UnitMatrix-method
##' imxCreateMatrix,ZeroMatrix-method
setGeneric("imxCreateMatrix", 
	function(.Object, labels, values, free, lbound, ubound, nrow, ncol, byrow, name, persist, ...) {
		return(standardGeneric("imxCreateMatrix"))
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
		lbound = "matrix", ubound = "matrix",
		.squareBrackets = "matrix", .persist = "logical",
		display = "character", dependencies = "integer",
	  "VIRTUAL"))

setMethod("imxCreateMatrix", "MxMatrix",
	function(.Object, labels, values, free, lbound, ubound, nrow, ncol, byrow, name, persist=TRUE, ...) {
		.Object <- populateMatrixSlot(.Object, "labels", labels, nrow, ncol)
		.Object <- populateMatrixSlot(.Object, "values", values, nrow, ncol)
		.Object <- populateMatrixSlot(.Object, "free", free, nrow, ncol)
		.Object <- populateMatrixSlot(.Object, "lbound", lbound, nrow, ncol)
		.Object <- populateMatrixSlot(.Object, "ubound", ubound, nrow, ncol)
		.Object@name <- name
    .Object@.persist <- persist
		return(.Object)
	}
)	

populateMatrixSlot <- function(object, slotName, vals, nr, nc) {
    lendat <- length(vals)
	if (lendat > 1 && (nr * nc) %% lendat != 0) {
		if (((lendat > nr) && (lendat %/% nr) * nr != lendat) ||
			((lendat < nr) && (nr %/% lendat) * lendat != nr))
				warning(paste("data length", lendat, "is not a sub-multiple",
					"or multiple of the number of rows", nr,
					"for argument", omxQuotes(slotName), "in",
					deparse(width.cutoff = 400L, imxLocateFunction("mxMatrix"))), call. = FALSE)
		else if (((lendat > nc) && (lendat %/% nc) * nc != lendat) ||
			     ((lendat < nc) && (nc %/% lendat) * lendat != nc))
				warning(paste("data length", lendat, "is not a sub-multiple",
					"or multiple of the number of columns", nc,
					"for argument", omxQuotes(slotName), "in", 
					deparse(width.cutoff = 400L, imxLocateFunction("mxMatrix"))), call. = FALSE)
	}
    dnSave <- dimnames(object)
	slot(object, slotName) <- suppressWarnings(matrix(vals, nr, nc))
    dimnames(object) <- dnSave
	return(object)
}

setMethod("imxVerifyMatrix", "MxMatrix",
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
		select <- .Object@.squareBrackets
		subs <- .Object@labels[select]
		lapply(subs, verifySquareBracket, .Object@name)
		if (imxSquareMatrix(.Object)) {
			verifySquare(.Object)
		}
	}
)

setMethod("imxSymmetricMatrix", "MxMatrix",
	function(.Object) { return(FALSE) }
)

setMethod("imxSquareMatrix", "MxMatrix",
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

setMethod("dim", "MxMatrix",
	function(x) {
		return(dim(x@values))
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
		labels <- as.matrix(x@labels[i, j, drop = drop])
		values <- as.matrix(x@values[i, j, drop = drop])		
		free <- as.matrix(x@free[i, j, drop = drop])
		lbound <- as.matrix(x@lbound[i, j, drop = drop])
		ubound <- as.matrix(x@ubound[i, j, drop = drop])
		type <- class(x)[[1]]
		nrow <- nrow(labels)
		ncol <- ncol(labels)
		dnames <- dimnames(values)
		newMatrix <- tryCatch(suppressWarnings(
			mxMatrix(type, nrow, ncol, free, values, labels, 
				lbound, ubound, FALSE, dnames, name)),
				error = function(e) mxMatrix("Full",
					nrow, ncol, free, values, labels,
					lbound, ubound, FALSE, dnames, name)) 
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
		if (!is.null(value)) {
			if (length(value) != 2) {
				msg <- paste("The 'dimnames' argument to MxMatrix",
					omxQuotes(x@name), "must have a length of 2")
				stop(msg, call. = FALSE)
			}
			if (!is.null(value[[1]]) && length(value[[1]]) != nrow(x) || 
				!is.null(value[[2]]) && length(value[[2]]) != ncol(x)) {
				msg <- paste("The MxMatrix object", omxQuotes(x@name), 
					"has specified dimnames with dimensions",
					length(value[[1]]), "x", length(value[[2]]), "but the matrix",
					"is of dimensions", nrow(x), "x", ncol(x))
				stop(msg, call. = FALSE)
			}
		}
		dimnames(x@values) <- value
		dimnames(x@free) <- value
		dimnames(x@labels) <- value
		dimnames(x@lbound) <- value
		dimnames(x@ubound) <- value
		return(x)
	}
)

setMethod("$", "MxMatrix", imxExtractSlot)

setReplaceMethod("$", "MxMatrix",
	function(x, name, value) {
        if(name %in% c("labels","values","free","lbound","ubound")) {
            return(populateMatrixSlot(x, name, value, nrow(x), ncol(x)))
        } 
		return(imxReplaceSlot(x, name, value, check=TRUE))
	}
)

setMethod("names", "MxMatrix", slotNames)

matrixTypes <- c("Diag", "Full", "Iden", "Lower", "Stand", "Sdiag", "Symm", "Unit", "Zero")
squareMatrices <- c("Diag", "Iden", "Lower", "Stand", "Sdiag", "Symm")

matrixCheckDims <- function(type, values, free, labels, lbound, ubound, nrow, ncol) {
	inputs <- list(values, free, labels, lbound, ubound)
	areMatrices <- sapply(inputs, is.matrix)
	theMatrices <- inputs[areMatrices] 
	if (length(theMatrices) > 1) {
		matches <- sapply(theMatrices, function(x) { identical(dim(x), dim(theMatrices[[1]])) })
		if (!all(matches)) {
			if (is.na(nrow) && is.na(ncol)) {
				stop(paste("Two or more matrix inputs have different dimensions.",
					"Use the 'nrow' and 'ncol' arguments in",
					deparse(width.cutoff = 400L, imxLocateFunction("mxMatrix"))), call. = FALSE)
			}
			if (!(type %in% squareMatrices) && (is.na(nrow) || is.na(ncol))) {
				stop(paste("Two or more matrix inputs have different dimensions.",
					"Use the 'nrow' and 'ncol' arguments in",
					deparse(width.cutoff = 400L, imxLocateFunction("mxMatrix"))), call. = FALSE)
			}
		}
	}
	if (is.na(nrow) && is.na(ncol)) {
		if(length(theMatrices) == 0) {
			stop(paste("You must specify 'nrow' and 'ncol' arguments in",
					deparse(width.cutoff = 400L, imxLocateFunction("mxMatrix"))), call. = FALSE)
		}
		nrow <- dim(theMatrices[[1]])[[1]]
		ncol <- dim(theMatrices[[1]])[[2]]
	}
	if ((is.na(nrow) || is.na(ncol)) && (type %in% squareMatrices)) {
		if (is.na(nrow)) nrow <- ncol
		if (is.na(ncol)) ncol <- nrow
	}
	return(c(nrow, ncol))
}

mxMatrix <- function(type = "Full", nrow = NA, ncol = NA, 
	free = FALSE, values = NA, labels = NA, 
	lbound = NA, ubound = NA, byrow = getOption('mxByrow'), 
	dimnames = NA, name = NA) {
	if (missing(dimnames) && !missing(values) && !is.null(dimnames(values))) {
		dimnames <- dimnames(values)
	}
	if (all.na(values)) { values <- as.numeric(values) }
	if (all.na(labels)) { labels <- as.character(labels) }
	if (all.na(lbound)) { lbound <- as.numeric(lbound) }
	if (all.na(ubound)) { ubound <- as.numeric(ubound) }
	if (all.na(nrow)) { nrow <- as.numeric(nrow) }
	if (all.na(ncol)) { ncol <- as.numeric(ncol) }
	if (single.na(name)) {
		name <- imxUntitledName()
	}
	imxVerifyName(name, 0)
	if (!is.character(name)) {
		stop(paste("'name' argument must",
			"be a character vector in", 
			deparse(width.cutoff = 400L, imxLocateFunction("mxMatrix"))), call. = FALSE)
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
				deparse(width.cutoff = 400L, imxLocateFunction("mxMatrix"))), call. = FALSE)
		} else if (is.na(nrow)) {
			nrow <- ncol
		} else if (is.na(ncol)) {
			ncol <- nrow
		}
	} else if (is.na(nrow) || is.na(ncol)) {
		stop(paste("Both 'nrow' and 'ncol'",
					"must be specified on a",
					"non-square matrix in",
					deparse(width.cutoff = 400L, imxLocateFunction("mxMatrix"))), call. = FALSE)
	}
	values <- as.numeric.preserve(values)
	lbound <- as.numeric.preserve(lbound)
	ubound <- as.numeric.preserve(ubound)
	typeName <- paste(type, "Matrix", sep="")
	newMatrix <- new(typeName)
	newMatrix <- imxCreateMatrix(newMatrix, labels, values, 
		free, lbound, ubound, nrow, ncol, byrow, name)
	if(length(dimnames) == 1 && is.na(dimnames)) {
	} else {
		dimnames(newMatrix) <- dimnames
	}
	imxVerifyMatrix(newMatrix)
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
		stop(paste(omxQuotes(name), "argument to mxMatrix function",
			"must be a scalar, a vector, or a matrix in",
			deparse(width.cutoff = 400L, imxLocateFunction("mxMatrix"))), call. = FALSE)
	}
}


matrixCheckErrors <- function(type, values, free, labels, lbound, ubound, nrow, ncol, name) {
	if (is.na(match(type, matrixTypes))) {
		stop(paste("Type must be one of:", 
			paste(matrixTypes, collapse=" "),
			"in", deparse(width.cutoff = 400L, imxLocateFunction("mxMatrix"))), call. = FALSE)
	}
	matrixCheckArgument(values, 'values')
	matrixCheckArgument(free, 'free')
	matrixCheckArgument(labels, 'labels')
	matrixCheckArgument(lbound, 'lbound')
	matrixCheckArgument(ubound, 'ubound')
	if (!is.numeric(values)) {
		stop(paste("'values' argument to mxMatrix function",
			"must be of numeric type in", 
			deparse(width.cutoff = 400L, imxLocateFunction("mxMatrix"))), call. = FALSE)
	}
	if (!is.logical(free)) {
		stop(paste("'free' argument to mxMatrix function",
			"must be of logical type in", 
			deparse(width.cutoff = 400L, imxLocateFunction("mxMatrix"))), call. = FALSE)
	}
	if (!is.character(labels)) {
		stop(paste("'labels' argument to mxMatrix function",
			"must be of character type in", 
			deparse(width.cutoff = 400L, imxLocateFunction("mxMatrix"))), call. = FALSE)
	}
	if (!is.numeric(lbound)) {
		stop(paste("'lbound' argument to mxMatrix function",
			"must be of numeric type in", 
			deparse(width.cutoff = 400L, imxLocateFunction("mxMatrix"))), call. = FALSE)
	}
	if (!is.numeric(ubound)) {
		stop(paste("'ubound' argument to mxMatrix function",
			"must be of numeric type in", 
			deparse(width.cutoff = 400L, imxLocateFunction("mxMatrix"))), call. = FALSE)
	}
	lapply(labels, imxVerifyReference, -2)
	if(any(is.na(free))) {
		stop(paste("'free' argument to mxMatrix function",
			"cannot contain an NA in",
			deparse(width.cutoff = 400L, imxLocateFunction("mxMatrix"))), call. = FALSE)
	}
	if (length(nrow) != 1 || !is.numeric(nrow)) {
		stop(paste("'nrow' argument to mxMatrix function",
			"must be either NA or a single numeric value in",
			deparse(width.cutoff = 400L, imxLocateFunction("mxMatrix"))), call. = FALSE)
	}
	if (length(ncol) != 1 || !is.numeric(ncol)) {
		stop(paste("'ncol' argument to mxMatrix function",
			"must be either NA or a single numeric value in",
			deparse(width.cutoff = 400L, imxLocateFunction("mxMatrix"))), call. = FALSE)
	}
}

generateParameterListHelper <- function(mxMatrix, result, matrixNumber, freeVarGroups) {

	free <- mxMatrix@free
	labels <- mxMatrix@labels
	lbound <- mxMatrix@lbound
	ubound <- mxMatrix@ubound
	isSymmetric <- imxSymmetricMatrix(mxMatrix)

	group <- c(0L)
	if (length(freeVarGroups)) for (gx in seq(2, length(freeVarGroups), 2)) {
		if (any(mxMatrix@name %in% freeVarGroups[[gx]])) {
			group <- union(group, freeVarGroups[[gx - 1L]])
		}
	}

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
			entry <- list(minBounds, maxBounds, group,
				      c(matrixNumber, row, col))
			if (isSymmetric && row != col) {
				entry[[length(entry) + 1]] <- c(matrixNumber, col, row)
			}
			result[[length(result)+1]] <- entry
			names(result)[[length(result)]] <- paste(mxMatrix@name, 
				"[", rows[i], ",", cols[i], "]", sep ="")
		} else if (length(grep(imxSeparatorChar, parameterName, fixed = TRUE)) == 0) {
			if (!is.null(result[[parameterName]])) {
				original <- result[[parameterName]]
				original[[3]] <- union(original[[3]], group)
				original[[length(original) + 1]] <- c(matrixNumber, row, col)
				if (isSymmetric && row != col) {
					original[[length(original) + 1]] <- c(matrixNumber, col, row)
				}
				result[[parameterName]] <- original
			} else {
				entry <- list(minBounds, maxBounds, group, c(matrixNumber, row, col))
				if (isSymmetric && row != col) {
					entry[[length(entry) + 1]] <- c(matrixNumber, col, row)
				}
				result[[parameterName]] <- entry
			}
		}
	}
	return(result)
}

matchDefinitionVariable <- function(parameterName) {
	# definition variables are of the form paste(modelname,"data",colname, sep=".")

	if (hasSquareBrackets(parameterName)) {
		components <- splitSubstitution(parameterName)
		return(grep(".data.", components[2:3], fixed=TRUE, value=TRUE))
	}

	return(grep(".data.", parameterName, fixed=TRUE, value=TRUE))
}

# Definition variables is a list:
# each entry of the list is a sublist of length 3 or more
# first entry of the sublist: data number
# second entry of the sublist: column number
# remaining entries of the sublist: c(matrix number, row, column)
matrixDefinitions <- function(flatModel, mxMatrix, result, matrixNumber) {
	labels <- mxMatrix@labels
	select <- !is.na(labels) & !mxMatrix@free
	if (all(select == FALSE)) {
		return(result)
	}
	parameterNames <- labels[select]
	rows <- row(labels)[select]
	cols <- col(labels)[select]
	for(px in 1:length(parameterNames)) {
		defVariables <- matchDefinitionVariable(parameterNames[px])
		if (length(defVariables) == 0) next

		row <- rows[px] - 1L
		col <- cols[px] - 1L
		for(i in 1:length(defVariables)) {
			defVariable <- defVariables[[i]]
			if (!is.null(result[[defVariable]])) {
				original <- result[[defVariable]]
				original[[length(original) + 1]] <- c(matrixNumber, row, col)
				result[[defVariable]] <- original
			} else {
				components <- unlist(strsplit(defVariable, imxSeparatorChar, fixed = TRUE))
				dataname <- paste(components[1:2], collapse=".")
				dataNumber <- match(dataname, names(flatModel@datasets))
				if (is.na(dataNumber)) {
					stop(paste("Dataset", omxQuotes(dataname), "not found.",
						   "Referred to by definition variable",
						   omxQuotes(defVariable), "in matrix", omxQuotes(mxMatrix@name)))
				}
				mxdata <- flatModel@datasets[[dataNumber]]
				if (!.hasSlot(mxdata, 'observed')) {
					stop(paste("Definition variable", omxQuotes(defVariable),
						   "in matrix", omxQuotes(mxMatrix@name),
						   "refers to", omxQuotes(dataname),
						   "but this mxData has no observed data"))
				}
				observed <- mxdata@observed
				columnNumber <- match(components[3], colnames(observed))
				if (is.na(columnNumber)) {
					stop(paste("Definition variable", omxQuotes(defVariable),
						   "in matrix", omxQuotes(mxMatrix@name),
						   "refers to column", omxQuotes(components[3]),
						   "but this column is not found"))
				}
				result[[defVariable]] <- list(dataNumber - 1L, columnNumber - 1L,
							      c(matrixNumber, row, col))
			}
		}
	}
	return(result)
}

generateMatrixValuesHelper <- function(mxMatrix) {
	return(mxMatrix@values)
}

displayMatrix <- function(mxMatrix) {
	type <- class(mxMatrix)[[1]]
	cat(type, omxQuotes(mxMatrix@name), '\n')
	cat("\n")
	nolabels <- all(is.na(mxMatrix@labels))
	if(nolabels == FALSE) {
		cat("$labels\n")
		print(mxMatrix@labels)
		cat("\n")
	} else {
		cat("$labels: No labels assigned.\n\n")
	}
	cat("$values\n")
	print(mxMatrix@values)
	cat("\n")
	noFree <- all(mxMatrix@free == FALSE)
	if(noFree == FALSE) {
		cat("$free\n")
		print(mxMatrix@free)
		cat("\n")
	} else {
		cat("$free: No free parameters.\n\n")
	}
	if(!all(is.na(mxMatrix@lbound))) {
		cat("$lbound\n")
		print(mxMatrix@lbound)
		cat("\n")
	} else {
		cat("$lbound: No lower bounds assigned.\n\n")
	}
	if(!all(is.na(mxMatrix@ubound))) {
		cat("$ubound\n")
		print(mxMatrix@ubound)
		cat("\n")
	} else {
		cat("$ubound: No upper bounds assigned.\n\n")
	}
}

setMethod("print", "MxMatrix", function(x,...) { displayMatrix(x) })
setMethod("show", "MxMatrix", function(object) { displayMatrix(object) })
