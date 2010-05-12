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

omxGetParameters <- function(model, indep = FALSE) {
	parameters <- lapply(model@matrices, getParametersHelper)
	names(parameters) <- NULL
	parameters <- unlist(parameters)
	if(indep) {
		submodels <- model@submodels
	} else {
		submodels <- omxDependentModels(model)
	}
	subparams <- lapply(submodels, omxGetParameters, indep)
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
	newlabels = NULL, lbound = NULL, ubound = NULL, indep = FALSE) {
	if (missing(labels) || !is.character(labels) || length(labels) == 0) {
		stop("'labels' argument must be a character vector")
	}
	if (any(is.na(labels))) {
		stop("'labels' argument must not contain NA values")	
	}
	if (any(duplicated(labels))) {
		stop("'labels' argument must not contain duplicate values")		
	}
	setParametersCheckVector(free, is.logical, 'free', 'logical')
	setParametersCheckVector(values, is.numeric, 'values', 'numeric')
	setParametersCheckVector(newlabels, is.character, 'newlabels', 'character')
	setParametersCheckVector(lbound, is.numeric, 'lbound', 'numeric')
	setParametersCheckVector(ubound, is.numeric, 'ubound', 'numeric')
	model@matrices <- lapply(model@matrices, setParametersHelper, 
		labels, free, values, newlabels, lbound, ubound)
	if(indep) {
		model@submodels <- lapply(model@submodels, omxSetParameters, 
			labels, free, values, newlabels, lbound, ubound, indep)
	} else {
		select <- omxDependentModels(model)
		select <- lapply(select, omxSetParameters, 
			labels, free, values, newlabels, lbound, ubound, indep)
		model@submodels <- c(select, omxIndependentModels(model))
	}
	return(model)
}

omxAssignFirstParameters <- function(model, indep = FALSE) {
	params <- omxGetParameters(model, indep)
	pnames <- names(params)
	model <- omxSetParameters(model, pnames[!is.na(pnames)], 
		values = params[!is.na(pnames)], indep = indep)
	return(model)
}

getParametersHelper <- function(amatrix) {
	free <- amatrix@free
	if (all(!free)) {
		return(numeric())
	}
	if (omxSymmetricMatrix(amatrix)) {
		triangle <- upper.tri(free, diag=TRUE)
		select <- free & triangle
	} else {
		select <- free
	}
	freeNames <- amatrix@labels[select]
	freeValues <- amatrix@values[select]
	names(freeValues) <- freeNames
	return(freeValues[!duplicated(freeNames, incomparables = NA)])
}


setParametersHelper <- function(amatrix, names, free, values, newlabels, lbound, ubound) {	
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


