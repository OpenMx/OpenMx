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

omxSetParameters <- function(model, names, values, indep = FALSE) {
	if(length(names) != length(values)) {
		stop("length of 'names' and 'values' must be identical")
	}
	model@matrices <- lapply(model@matrices, setParametersHelper, names, values)
	if(indep) {
		model@submodels <- lapply(model@submodels, omxSetParameters, names, values, indep)
	} else {
		select <- omxDependentModels(model)
		select <- lapply(select, omxSetParameters, names, values, indep)
		model@submodels <- c(select, omxIndependentModels(model))
	}
	return(model)
}

omxAssignFirstParameters <- function(model, indep = FALSE) {
	params <- omxGetParameters(model, indep)
	pnames <- names(params)
	model <- omxSetParameters(model, pnames[!is.na(pnames)], params[!is.na(pnames)], indep)
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


setParametersHelper <- function(amatrix, names, values) {	
	labels <- amatrix@labels
	locations <- which(labels %in% names)
	indices <- match(labels[locations], names)
	amatrix@values[locations] <- values[indices]
	return(amatrix)
}


