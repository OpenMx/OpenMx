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

shareData <- function(model) {
	return(shareDataHelper(model, model@data))
}

# If we reach an independent model and this model
# does not have any data, then give the current "default"
# data to this model. As we move down the tree, update
# the current "defult" data.
shareDataHelper <- function(model, default) {
	if((is.null(model@data)) && (model@independent == TRUE)) {
		model@data <- default
	} else if (!is.null(model@data)) {
		default <- model@data
	}
	submodels <- lapply(model@submodels, 
		shareDataHelper, default)
	model@submodels <- submodels
	return(model)
}

getAllIndependents <- function(model) {
    submodels <- model@submodels
    if(length(submodels) == 0) return(submodels)
    select <- sapply(submodels, function(x) { x@independent })
	retval <- submodels[select]
	recursive <- lapply(submodels[!select], getAllIndependents)
	recursive <- unlist(recursive, recursive = TRUE)
	retval <- c(retval, recursive)
    return(retval)
}

freezeMatrix <- function(mxMatrix) {
	rows <- nrow(mxMatrix)
	cols <- ncol(mxMatrix)
	mxMatrix@free <- matrix(FALSE, rows, cols)	
	return(mxMatrix)
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
		if (length(objective@result) > 0) {
			newMatrix <- mxMatrix(values = objective@result)
			newMatrix@name <- objective@name
			model[[objective@name]] <- newMatrix
		}
	}
	return(model)
}

imxFreezeModel <- function(model) {
	model <- freezeObjective(model)
	model@matrices <- lapply(model@matrices, freezeMatrix)
	algebras <- lapply(model@algebras, freezeAlgebra)
	algebras <- algebras[!sapply(algebras, is.null)]
	model@matrices <- append(model@matrices, algebras)
	model@algebras <- list()
	model@constraints <- list()
	model@submodels <- lapply(model@submodels, imxFreezeModel)
	return(model)
}

imxFlattenModel <- function(model, namespace) {
	flatModel <- new("MxFlatModel", model, list(), list())
	name <- model@name
	model@objective <- namespaceConvertObjective(model@objective, name, namespace)
	model@data <- namespaceConvertData(model@data, name)
	flatModel@objective <- model@objective
	defaultData <- model@data
	flatModel@data <- defaultData
	if (is.null(defaultData)) {
		defaultDataName <- NULL
	} else {
		defaultDataName <- defaultData@name
	}
	flatModel@matrices <- lapply(model@matrices, namespaceConvertMatrix, name, defaultDataName, namespace)
	names(flatModel@matrices) <- imxExtractNames(flatModel@matrices)	
	flatModel@algebras <- collectComponents(model, namespace, "algebras", namespaceConvertAlgebra)
	flatModel@constraints <- collectComponents(model, namespace, "constraints", namespaceConvertConstraint)	
	flatModel@intervals <- collectComponents(model, namespace, "intervals", namespaceConvertInterval)		
	flatModel <- flattenModelHelper(model, flatModel, defaultData, namespace)
	flatModel@submodels <- list()
	return(flatModel)
}

collectComponents <- function(model, namespace, slotName, convertFunction) {
	components <- collectComponentsHelper(model, namespace, slotName, convertFunction)
	if (slotName != "intervals") {
		names(components) <- imxExtractNames(components)
	}
	return(components)
}

collectComponentsHelper <- function(model, namespace, slotName, convertFunction) {
	components <- lapply(slot(model, slotName), convertFunction, model@name, namespace)
	if (length(model@submodels) > 0) {
		submodel_components <- lapply(model@submodels, collectComponents, namespace, slotName, convertFunction)
		submodel_components <- unlist(submodel_components, recursive = FALSE, use.names = FALSE)
		components <- append(components, submodel_components)
	}
	return(components)
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
			if (is.null(defaultData)) {
				defaultDataName <- NULL
			} else {
				defaultDataName <- defaultData@name
			}
			if (is.null(submodel@data)) {
				submodel@matrices <- lapply(submodel@matrices, namespaceConvertMatrix,
					name, defaultDataName, namespace)
			} else {
				submodel@matrices <- lapply(submodel@matrices, namespaceConvertMatrix,
					name, submodel@data@name, namespace)
			}
			names(submodel@matrices) <- imxExtractNames(submodel@matrices)
			flatModel@matrices    <- append(flatModel@matrices, submodel@matrices)
			flatModel <- flattenModelHelper(submodel, flatModel, defaultData, namespace)
		}
	}
	return(flatModel)
}

imxDependentModels <- function(model) {
        retval <- model@submodels
        if(length(retval) == 0) return(retval)
        retval <- retval[sapply(retval, function(x) { !x@independent })]
        return(retval)
}

imxIndependentModels <- function(model) {
        retval <- model@submodels
        if(length(retval) == 0) return(retval)
        retval <- retval[sapply(retval, function(x) { x@independent })]
        return(retval)
}


imxReplaceModels <- function(model, replacements) {
	if (length(replacements) == 0) return(model)
	names(replacements) <- imxExtractNames(replacements)
	return(replaceModelsHelper(model, replacements))
}

replaceSubmodels <- function(target, replacements) {
	retval <- replacements[[target@name]]
	if (is.null(retval)) {
		return(target)
	} else {
		return(retval)
	}
}

replaceModelsHelper <- function(model, replacements) {
	submodels <- model@submodels
	if (length(submodels) == 0) return(model)
	submodels <- lapply(submodels, replaceSubmodels, replacements)
	submodels <- lapply(submodels, replaceModelsHelper, replacements)
	model@submodels <- submodels
	return(model)
}
