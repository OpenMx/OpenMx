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

omxGetParameters <- function(model, indep = FALSE, free = c(TRUE, FALSE, NA)) {
	if (identical(free, c(TRUE, FALSE, NA))) {
		free <- TRUE
	}
	if (!is.logical(free) || length(free) != 1) {
		stop("argument 'free' must be a 'TRUE', 'FALSE', or NA")
	}
	parameters <- lapply(model@matrices, getParametersHelper, free)
	names(parameters) <- NULL
	parameters <- unlist(parameters)
	if(indep) {
		submodels <- model@submodels
	} else {
		submodels <- imxDependentModels(model)
	}
	subparams <- lapply(submodels, omxGetParameters, indep, free)
	names(subparams) <- NULL
	subparams <- unlist(subparams)
	parameters <- c(parameters, subparams)
	parameters <- parameters[!duplicated(names(parameters), incomparables = NA)]
	return(parameters)
}

setParametersCheckVector <- function(values, test, argname, typename) {
	if (is.null(values)) return()
	if (!test(values)) {
		stop(paste(omxQuotes(argname), 
			"argument must either be NA or a",
			typename, "vector"), call. = FALSE)
	}
}


omxSetParameters <- function(model, labels, free = NULL, values = NULL,
	newlabels = NULL, lbound = NULL, ubound = NULL, indep = FALSE,
	strict = TRUE) {
	if (missing(labels) || !is.character(labels) || length(labels) == 0) {
		stop("'labels' argument must be a character vector")
	}
	if (any(is.na(labels))) {
		stop("'labels' argument must not contain NA values")	
	}
	if (any(duplicated(labels))) {
		stop("'labels' argument must not contain duplicate values")		
	}
	if (strict) {
		pnames <- names(omxGetParameters(model, indep, NA))
		missing <- setdiff(labels, pnames)
		if (length(missing) > 0) {
			msg <- paste("The following labels are",
				"not present in the model",
				"(use 'strict' = FALSE to ignore):",
				omxQuotes(missing))
			stop(msg)
		}
	}
	if (is.vector(lbound) && length(lbound) > 0 && all(sapply(lbound, is.na))) {
		lbound <- as.numeric(lbound)
	}
	if (is.vector(ubound) && length(ubound) > 0 && all(sapply(ubound, is.na))) {
		ubound <- as.numeric(ubound)
	}
	setParametersCheckVector(free, is.logical, 'free', 'logical')
	setParametersCheckVector(values, is.numeric, 'values', 'numeric')
	setParametersCheckVector(newlabels, is.character, 'newlabels', 'character')
	setParametersCheckVector(lbound, is.numeric, 'lbound', 'numeric')
	setParametersCheckVector(ubound, is.numeric, 'ubound', 'numeric')
	return(setParametersHelper(model, labels, free, values,
		newlabels, lbound, ubound, indep))
}

setParametersHelper <- function(model, labels, free, values,
	newlabels, lbound, ubound, indep) {
	model@matrices <- lapply(model@matrices, setParametersMatrix, 
		labels, free, values, newlabels, lbound, ubound)
	if(indep) {
		if (length(model@submodels) == 0) {
			return(model)
		}
		model@submodels <- lapply(model@submodels, setParametersHelper, 
			labels, free, values, newlabels, lbound, ubound, indep)
	} else {
		select <- imxDependentModels(model)
		if (length(select) == 0) {
			return(model)
		}
		select <- lapply(select, setParametersHelper, 
			labels, free, values, newlabels, lbound, ubound, indep)
		model@submodels <- c(select, imxIndependentModels(model))
	}
	return(model)
}

extractFirst <- function(x) {
	return(x[[1]])
}

extractSecond <- function(x) {
	return(x[[2]])
}

omxNameAnonymousParameters <- function(model, indep = FALSE) {
	rows <- lapply(model@matrices, getAnonymousRows)
    cols <- lapply(model@matrices, getAnonymousCols)
    newnames <- mapply(getAnonymousNames, rows)
	model@matrices <- mapply(assignAnonymousNames, model@matrices, rows, cols, newnames)
	newnames <- unlist(newnames)
	if(indep) {
		if (length(model@submodels) == 0) {
			return(list(model, newnames))
		}
		pairs <- lapply(model@submodels, omxNameAnonymousParameters, indep)
		submodels <- lapply(pairs, extractFirst)
		subnames <- unlist(lapply(pairs, extractSecond))
		names(submodels) <- names(model@submodels)
		model@submodels <- submodels
		newnames <- c(newnames, subnames)
	} else {
		select <- imxDependentModels(model)
		if (length(select) == 0) {
			return(list(model, newnames))
		}
		pairs <- lapply(select, omxNameAnonymousParameters, indep)
        submodels <- lapply(pairs, extractFirst)
        subnames <- unlist(lapply(pairs, extractSecond))
        names(submodels) <- names(select)
		model@submodels <- c(select, imxIndependentModels(model))
        newnames <- c(newnames, subnames)
	}
	return(list(model, newnames))
}

getAnonymousRows <- function(matrix) {
	select <- matrix@free & is.na(matrix@labels)
	if (imxSymmetricMatrix(matrix)) {
		select <- select & upper.tri(matrix@labels, diag=TRUE)
	}
	return(row(matrix@free)[select])
}

getAnonymousCols <- function(matrix) {
	select <- matrix@free & is.na(matrix@labels)
	if (imxSymmetricMatrix(matrix)) {
		select <- select & upper.tri(matrix@labels, diag=TRUE)
	}
	return(col(matrix@free)[select])
}

assignAnonymousNames <- function(matrix, rows, cols, newnames) {
	symmetry <- imxSymmetricMatrix(matrix)
	if (length(rows) > 0) {
		for(i in 1:length(rows)) {
			row <- rows[[i]]
			col <- cols[[i]]
			newname <- newnames[[i]]
			matrix@labels[row,col] <- newname
			if (symmetry) {
				matrix@labels[col,row] <- newname
			}
		}
	}
	return(matrix)
}

getAnonymousNames <- function(rows) {
	return(replicate(length(rows), imxUntitledName()))
}

omxAssignFirstParameters <- function(model, indep = FALSE) {
	params <- omxGetParameters(model, indep)
	pnames <- names(params)
	model <- omxSetParameters(model, pnames[!is.na(pnames)], 
		values = params[!is.na(pnames)], indep = indep)
	return(model)
}

getParametersHelper <- function(amatrix, selection) {
	if (single.na(selection)) {
		select <- amatrix@free | !apply(amatrix@labels, c(1,2), is.na)
	} else if (selection) {
		select <- amatrix@free
	} else {
		select <- !amatrix@free & !apply(amatrix@labels, c(1,2), is.na)
	}
	if (all(!select)) {
		return(numeric())
	}
	if (imxSymmetricMatrix(amatrix)) {
		triangle <- upper.tri(select, diag=TRUE)
		select <- select & triangle
	} 
	theNames <- amatrix@labels[select]
	theValues <- amatrix@values[select]
	names(theValues) <- theNames
	return(theValues[!duplicated(theNames, incomparables = NA)])
}


setParametersMatrix <- function(amatrix, names, free, values, newlabels, lbound, ubound) {	
	labels <- amatrix@labels
	locations <- which(labels %in% names)
	indices <- match(labels[locations], names)
	if (!is.null(free)) {
		index2 <- ((indices - 1) %% length(free)) + 1
		amatrix@free[locations] <- as.logical(free[index2])
	}
	if (!is.null(values)) {
		index2 <- ((indices - 1) %% length(values)) + 1
		amatrix@values[locations] <- as.numeric(values[index2])
	}
	if (!is.null(newlabels)) {
		index2 <- ((indices - 1) %% length(newlabels)) + 1
		amatrix@labels[locations] <- as.character(newlabels[index2])
	}
	if (!is.null(lbound)) {
		index2 <- ((indices - 1) %% length(lbound)) + 1
		amatrix@lbound[locations] <- as.numeric(lbound[index2])
	}
	if (!is.null(ubound)) {
		index2 <- ((indices - 1) %% length(ubound)) + 1
		amatrix@ubound[locations] <- as.numeric(ubound[index2])
	}
	return(amatrix)
}


