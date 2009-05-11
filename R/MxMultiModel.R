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

shareData <- function(model) {
	return(shareDataHelper(model, model@data))
}

# If we reach an independent model and this model
# does not have any data, then give the current "default"
# data to this model. As we move down the true, update
# the current "defult" data.
shareDataHelper <- function(model, default) {
	if((is.null(model@data)) && (model@independent == TRUE)) {
		model@data <- default
	} else if (!is.null(model@data)) {
		default <- model@data
	}
	submodels <- lapply(model@submodels, function(x)
		{ shareDataHelper(x, default) })
	model@submodels <- submodels
	return(model)
}

omxGetIndependents <- function(model) {
	return(getIndependentsHelper(model, list()))
}

getIndependentsHelper <- function(model, lst) {
	indep <- lapply(model@submodels, function(x) {
		if(x@independent == TRUE) {x}
		else {NULL}
	})
	res <- append(lst, indep[!sapply(indep, is.null)])
	if (length(model@submodels) > 0) {
		for(i in 1:length(model@submodels)) {
			res <- getIndependentsHelper(model@submodels[[i]], res)
		}
	}
	return(res)
}

freezeMatrix <- function(mxMatrix) {
	type <- class(mxMatrix)[[1]]
	rows <- nrow(mxMatrix)
	cols <- ncol(mxMatrix)	
	return(new(type, mxMatrix@name, mxMatrix@values, 
		FALSE, mxMatrix@labels, rows, cols))
}

freezeAlgebra <- function(mxAlgebra) {
	if(is.null(mxAlgebra@result)) return(NULL)
	res <- mxMatrix(values = mxAlgebra@result, 
		name = mxAlgebra@name)
	return(res)
}

freezeObjective <- function(model) {
	objective <- model@objective
	if (!is.null(objective)) {
		model[[objective@name]] <- NULL
		if (!is.null(objective@result)) {
			model[[objective@name]] <- mxMatrix(
				values = objective@result,
				name = objective@name)
		}
	}
	return(model)
}

omxFreezeModel <- function(model) {
	model <- freezeObjective(model)
	model@matrices <- lapply(model@matrices, freezeMatrix)
	algebras <- lapply(model@algebras, freezeAlgebra)
	algebras <- algebras[!sapply(algebras, is.null)]
	model@matrices <- append(model@matrices, algebras)
	model@algebras <- list()
	model@submodels <- lapply(model@submodels, omxFreezeModel)
	return(model)
}

omxFlattenModel <- function(model, namespace) {
	flatModel <- new("MxFlatModel", model, list(), list())
	name <- model@name
	model@objective <- namespaceConvertObjective(model@objective, name, namespace)
	model@data <- namespaceConvertData(model@data, name)
	flatModel@objective <- model@objective
	defaultData <- model@data
	flatModel@data <- defaultData
	flatModel@matrices <- lapply(model@matrices, 
		function(x) { namespaceConvertMatrix(x, name, defaultData@name, namespace) })
	flatModel@algebras <- lapply(model@algebras, 
		function(x) { namespaceConvertAlgebra(x, name, namespace) })
	flatModel@constraints <- lapply(model@constraints,
		function(x) { namespaceConvertConstraint(x, name, namespace) })
	names(flatModel@matrices) <- entityExtractNames(flatModel@matrices)
	names(flatModel@algebras) <- entityExtractNames(flatModel@algebras)
	names(flatModel@constraints) <- entityExtractNames(flatModel@constraints)
	flatModel <- flattenModelHelper(model, flatModel, defaultData, namespace)
	flatModel@submodels <- list()
	return(flatModel)
}

flattenModelHelper <- function(model, flatModel, defaultData, namespace) {
	if (!is.null(model@objective)) {
		if(is.na(model@objective@data) && is.null(model@data) && !is.null(defaultData)) {
			model@objective@data <- defaultData@name
		} else if (is.na(model@objective@data) && !is.null(model@data)) {
			model@objective@data <- model@data@name
		}
		flatModel@objectives[[model@objective@name]] <- model@objective
	}
	if (!is.null(model@data)) {
		flatModel@datasets[[model@data@name]] <- model@data
	}
	if (is.null(defaultData)) {
		defaultData <- model@data
	}
	if (length(model@submodels) > 0) {
		for(i in 1:length(model@submodels)) {
			submodel <- model@submodels[[i]]
			name <- submodel@name
			submodel@data <- namespaceConvertData(submodel@data, name)
			submodel@objective <- namespaceConvertObjective(submodel@objective, name, namespace)
			if (is.null(submodel@data)) {
				submodel@matrices <- lapply(submodel@matrices, 
					function(x) { 
						namespaceConvertMatrix(x, name, defaultData@name, namespace) })
			} else {
				submodel@matrices <- lapply(submodel@matrices, 
					function(x) { 
						namespaceConvertMatrix(x, name, submodel@data@name, namespace) })
			}
			submodel@algebras <- lapply(submodel@algebras, 
				function(x) { namespaceConvertAlgebra(x, name, namespace) })
			submodel@constraints <- lapply(submodel@constraints, 
				function(x) { namespaceConvertConstraint(x, name, namespace) })
			names(submodel@matrices) <- entityExtractNames(submodel@matrices)
			names(submodel@algebras) <- entityExtractNames(submodel@algebras)
			names(submodel@constraints) <- entityExtractNames(submodel@constraints)
			flatModel@matrices    <- append(flatModel@matrices, submodel@matrices)
			flatModel@algebras    <- append(flatModel@algebras, submodel@algebras)
			flatModel@constraints <- append(flatModel@constraints, 
				submodel@constraints) 
			flatModel <- flattenModelHelper(submodel, flatModel, defaultData)
		}
	}
	return(flatModel)
}

omxReplaceModels <- function(model, replacements) {
	mnames <- names(model@submodels)
	rnames <- names(replacements)
	inames <- intersect(mnames, rnames)
	if (length(inames) > 0) {
		for(i in 1:length(inames)) {
			name <- inames[[i]]
			model[[name]] <- replacements[[name]]
		}
	}
	if (length(model@submodels) > 0) {
		for(i in 1:length(model@submodels)) {
			model@submodels[[i]] <- 
				omxReplaceModels(model@submodels[[i]], replacements)
		}
	}
	return(model)
}
