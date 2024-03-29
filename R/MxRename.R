#
#   Copyright 2007-2021 by the individuals mentioned in the source code history
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

getAllModelNames <- function(model) {
	ret <- c(model@name)
	if (length(model@submodels)) ret <- c(ret, sapply(model@submodels, getAllModelNames))
	ret
}

mxRename <- function(model, newname, oldname = NA) {
  warnModelCreatedByOldVersion(model)
	if( !is.character(newname)) {
		stop("'newname' argument is not a character string")
	}
	if( !(is.na(oldname) || is.character(oldname))) {
		stop("'oldname' argument is not either NA or character string")
	}
	imxVerifyName(newname, 0)
	if (is.na(oldname)) {
		oldname <- model@name
	}
	if (newname == oldname) return(model)
	existing <- getAllModelNames(model)
	if (newname %in% existing) {
		stop(paste("There is already a model named", omxQuotes(newname)))
	}
	model <- propagateModelName(model, oldname, newname)
	return(model)
}

propagateModelName <- function(model, oldname, newname) {
	if(model@name == oldname) {
		model@name <- newname
	}
	model@matrices <- lapply(model@matrices, renameMatrix, oldname, newname)
	model@algebras <- lapply(model@algebras, renameAlgebra, oldname, newname)
	model@constraints <- lapply(model@constraints, renameConstraint, oldname, newname)
	model@intervals <- lapply(model@intervals, renameConfidenceIntervals, oldname, newname)
	model@expectation <- genericExpRename(model@expectation, oldname, newname)
	model@fitfunction <- genericFitRename(model@fitfunction, oldname, newname)
	model@submodels <- lapply(model@submodels, propagateModelName, oldname, newname)
	model@output <- renameModelOutput(model@output, oldname, newname)
	names(model@intervals) <- imxExtractReferences(model@intervals)
	names(model@submodels) <- imxExtractNames(model@submodels)
	return(model)
}

renameReference <- function(reference, oldname, newname) {
  if (length(reference) == 0) return(reference)
	if (is.na(reference)) {
		return(reference)
	}
	components <- unlist(strsplit(reference, imxSeparatorChar, fixed = TRUE))
	if (length(components) == 2 && components[[1]] == oldname) {
		return(paste(newname, components[[2]], sep = imxSeparatorChar))
  } else if (length(components) == 1 && components[[1]] == oldname) {
    # Don't have enough context to know whether this is a model name or not.
		reference
	} else {
		reference
	}
}

renameUnqualifiedReference <- function(reference, oldname, newname) {
  if (length(reference) == 0) return(reference)
	if (is.na(reference)) {
		return(reference)
	}
	components <- unlist(strsplit(reference, imxSeparatorChar, fixed = TRUE))
	if (length(components) == 2 && components[[1]] == oldname) {
		return(paste(newname, components[[2]], sep = imxSeparatorChar))
  } else if (length(components) == 1 && components[[1]] == oldname) {
    newname  # assume it's a model name
	} else {
		reference
	}
}

renameModelOutput <- function(output, oldname, newname) {
	if(is.null(output)) {
		return(output)
	}
	if(!is.null(output$confidenceIntervals)) {
		names <- dimnames(output$confidenceIntervals)
		rownames <- lapply(names[[1]], renameReference, oldname, newname)
		dimnames(output$confidenceIntervals) <- list(rownames, names[[2]])
	}
	if(!is.null(output$confidenceIntervalCodes)) {
		names <- dimnames(output$confidenceIntervalCodes)
		rownames <- lapply(names[[1]], renameReference, oldname, newname)
		dimnames(output$confidenceIntervalCodes) <- list(rownames, names[[2]])
	}
	return(output)
}

renameMatrix <- function(matrix, oldname, newname) {
	matrix@labels <- apply(matrix@labels, c(1,2), renameReference, oldname, newname)
	return(matrix)
}

renameConstraint <- function(constraint, oldname, newname) {
	constraint@name <- renameReference(constraint@name, oldname, newname)
	constraint@formula <- renameFormula(constraint@formula, oldname, newname)
	constraint@alg1 <- renameReference(constraint@alg1, oldname, newname)
	constraint@alg2 <- renameReference(constraint@alg2, oldname, newname)
	return(constraint)
}


renameSymbol <- function(symbol, oldname, newname) {
	if (is.numeric(symbol) || is.character(symbol)) {
		return(symbol)
	} else {
		return(as.symbol(renameReference(as.character(symbol), oldname, newname)))
	}
}

renameAlgebra <- function(algebra, oldname, newname) {
	algebra@formula <- renameFormula(algebra@formula, oldname, newname)
	return(algebra)
}

renameConfidenceIntervals <- function(interval, oldname, newname) {
	interval@reference <- renameReference(interval@reference, oldname, newname)
	return(interval)
}

renameFormula <- function(formula, oldname, newname) {
	len <- length(formula)
	if (len == 0) {
		stop("mxRename has reached an invalid state")
	} else if (len == 1) {
		formula <- renameSymbol(formula, oldname, newname)
	} else {
		tail <- formula[-1]
		select <- sapply(tail, function(x) { !identical(as.character(x), '') })
		tail[select] <- lapply(tail[select], renameFormula, oldname, newname)
		formula[-1] <- tail
	}
	return(formula)
}
