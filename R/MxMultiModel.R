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
	flatModel@objective <- namespaceConvertObjective(model@objective, name, character(), namespace)
	defaultData <- namespaceConvertData(model@data, name)
	flatModel@data <- defaultData
	flatModel@matrices <- collectMatrices(model, namespace, character(), 
		defaultData)
	flatModel@algebras <- collectComponents(model, namespace, "algebras", 
		character(), namespaceConvertAlgebra)
	flatModel@constraints <- collectComponents(model, namespace, "constraints", 
		character(), namespaceConvertConstraint)	
	flatModel@intervals <- collectComponents(model, namespace, "intervals", 
		character(), namespaceConvertInterval)
	flatModel@datasets <- collectDatasets(model)
	flatModel@objectives <- collectObjectives(model, namespace, character(), 
		defaultData)
	flatModel@submodels <- list()
	return(flatModel)
}

updateContext <- function(context, modelname) {
	if (length(context) == 0) {
		return(modelname)
	} else {
		return(paste(context, modelname, sep = '.'))
	}
}

collectComponents <- function(model, namespace, slotName, context, convertFunction) {
	components <- collectComponentsHelper(model, namespace, slotName, context, convertFunction)
	if (slotName != "intervals") {
		names(components) <- imxExtractNames(components)
	}
	return(components)
}

collectComponentsHelper <- function(model, namespace, slotName, context, convertFunction) {
	components <- lapply(slot(model, slotName), convertFunction, model@name, context, namespace)
	if (length(model@submodels) > 0) {
		context <- updateContext(context, model@name)
		submodel_components <- lapply(model@submodels, collectComponentsHelper, 
			namespace, slotName, context, convertFunction)
		submodel_components <- unlist(submodel_components, recursive = FALSE, use.names = FALSE)
		components <- append(components, submodel_components)
	}
	return(components)
}

collectDatasets <- function(model) {
	datasets <- collectDatasetsHelper(model)
	names(datasets) <- imxExtractNames(datasets)
	return(datasets)
}


collectDatasetsHelper <- function(model) {
	modeldata <- model@data
	if (!is.null(modeldata)) {
		modeldata <- namespaceConvertData(modeldata, model@name)
		retval <- list(modeldata)
	} else {
		retval <- list()
	}
	if (length(model@submodels) > 0) {
		submodel_datasets <- lapply(model@submodels, collectDatasetsHelper)		
		submodel_datasets <- unlist(submodel_datasets, recursive = FALSE, use.names = FALSE)
		retval <- append(retval, submodel_datasets)
	}
	return(retval)
}

collectMatrices <- function(model, namespace, context, defaultData) {
	matrices <- collectMatricesHelper(model, namespace, context, defaultData)
	names(matrices) <- imxExtractNames(matrices)
	return(matrices)
}

collectMatricesHelper <- function(model, namespace, context, defaultData) {
	modeldata <- namespaceConvertData(model@data, model@name)
	if (is.null(defaultData)) {
		defaultData <- modeldata
	} 
	if (is.null(defaultData)) {
		defaultDataName <- NULL
	} else {
		defaultDataName <- defaultData@name
	}
	if (is.null(modeldata)) {
		retval <- lapply(model@matrices, namespaceConvertMatrix,
			model@name, defaultDataName, context, namespace)
	} else {
		retval <- lapply(model@matrices, namespaceConvertMatrix,
			model@name, modeldata@name, context, namespace)
	}
	if (length(model@submodels) > 0) {
		context <- updateContext(context, model@name)
		submodel_matrices <- lapply(model@submodels, collectMatricesHelper, namespace, context, defaultData)		
		submodel_matrices <- unlist(submodel_matrices, recursive = FALSE, use.names = FALSE)
		retval <- append(retval, submodel_matrices)
	}
	return(retval)
}

collectObjectives <- function(model, namespace, context, defaultData) {
	objectives <- collectObjectivesHelper(model, namespace, context, defaultData)
	names(objectives) <- imxExtractNames(objectives)
	return(objectives)
}


collectObjectivesHelper <- function(model, namespace, context, defaultData) {
	objective <- namespaceConvertObjective(model@objective, model@name, context, namespace)
	modeldata <- namespaceConvertData(model@data, model@name)	
	if (is.null(defaultData)) {
		defaultData <- modeldata
	} 	
	if (!is.null(objective)) {
		if(is.na(objective@data) && is.null(modeldata) && !is.null(defaultData)) {
			objective@data <- defaultData@name
		} else if (is.na(objective@data) && !is.null(modeldata)) {
			objective@data <- modeldata@name
		}
		retval <- list(objective)
	} else {
		retval <- list()
	}
	if (length(model@submodels) > 0) {
		context <- updateContext(context, model@name)
		submodel_objectives <- lapply(model@submodels, collectObjectivesHelper, namespace, context, defaultData)		
		submodel_objectives <- unlist(submodel_objectives, recursive = FALSE, use.names = FALSE)
		retval <- append(retval, submodel_objectives)
	}
	return(retval)	
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
