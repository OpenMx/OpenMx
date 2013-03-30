#
#   Copyright 2007-2013 The OpenMx Project
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

freezeFitfunction <- function(model) {
	fitfunction <- model@fitfunction
	if (!is.null(fitfunction)) {
		model[[fitfunction@name]] <- NULL
		if (length(fitfunction@result) > 0) {
			newMatrix <- mxMatrix(values = fitfunction@result)
			newMatrix@name <- fitfunction@name
			model[[fitfunction@name]] <- newMatrix
		}
	}
	return(model)
}

imxFreezeModel <- function(model) {
	model <- freezeFitfunction(model)
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
	flatModel <- new("MxFlatModel", model, list(), list(), list())
	name <- model@name
	flatModel@fitfunction <- namespaceConvertFitFunction(model@fitfunction, name, namespace)
	flatModel@expectation <- namespaceConvertExpectation(model@expectation, name, namespace)
	defaultData <- namespaceConvertData(model@data, name)
	flatModel@data <- defaultData
	flatModel@matrices <- collectMatrices(model, namespace, defaultData)
	flatModel@algebras <- collectComponents(model, namespace, "algebras", namespaceConvertAlgebra)
	flatModel@constraints <- collectComponents(model, namespace, "constraints", namespaceConvertConstraint)	
	flatModel@intervals <- collectComponents(model, namespace, "intervals", namespaceConvertInterval)
	flatModel@datasets <- collectDatasets(model)
	flatModel@fitfunctions <- collectFitFunctions(model, namespace, defaultData)
	flatModel@expectations <- collectExpectations(model, namespace, defaultData)
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
		submodel_components <- lapply(model@submodels, collectComponentsHelper, namespace, slotName, convertFunction)
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

collectMatrices <- function(model, namespace, defaultData) {
	matrices <- collectMatricesHelper(model, namespace, defaultData)
	names(matrices) <- imxExtractNames(matrices)
	return(matrices)
}

collectMatricesHelper <- function(model, namespace, defaultData) {
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
			model@name, defaultDataName, namespace)
	} else {
		retval <- lapply(model@matrices, namespaceConvertMatrix,
			model@name, modeldata@name, namespace)
	}
	if (length(model@submodels) > 0) {
		submodel_matrices <- lapply(model@submodels, collectMatricesHelper, namespace, defaultData)		
		submodel_matrices <- unlist(submodel_matrices, recursive = FALSE, use.names = FALSE)
		retval <- append(retval, submodel_matrices)
	}
	return(retval)
}

collectFitFunctions <- function(model, namespace, defaultData) {
	fitfunctions <- collectFitFunctionsHelper(model, namespace, defaultData)
	names(fitfunctions) <- imxExtractNames(fitfunctions)
	return(fitfunctions)
}

collectExpectations <- function(model, namespace, defaultData) {
	expectations <- collectExpectationsHelper(model, namespace, defaultData)
	if (length(expectations) == 0) return(list())
	names(expectations) <- imxExtractNames(expectations)
	return(expectations)
}

collectExpectationsHelper <- function(model, namespace, defaultData) {
	expectation <- namespaceConvertExpectation(model@expectation, model@name, namespace)
	modeldata <- namespaceConvertData(model@data, model@name)	
	if (is.null(defaultData)) {
		defaultData <- modeldata
	} 	
	container <- character(0)
	if (!is.null(expectation)) {
		if(is.na(expectation@data) && is.null(modeldata) && !is.null(defaultData)) {
			expectation@data <- defaultData@name
		} else if (is.na(expectation@data) && !is.null(modeldata)) {
			expectation@data <- modeldata@name
		}
		container <- expectation@name
	}
	submodel_expectations <- c()
	if (length(model@submodels) > 0) {
		submodel_expectations <- lapply(model@submodels, collectExpectationsHelper, namespace, defaultData)		
		submodel_expectations <- unlist(submodel_expectations, recursive = FALSE, use.names = FALSE)
		submodel_expectations <- lapply(submodel_expectations, function (e) {
			e@container <- container
			e
		})
		if (!is.null(expectation)) {
			expectation@submodels <- unlist(sapply(submodel_expectations, function(e) e@name))
		}
	}
	return(c(expectation, submodel_expectations))
}

collectFitFunctionsHelper <- function(model, namespace, defaultData) {
	fitfunction <- namespaceConvertFitFunction(model@fitfunction, model@name, namespace)
	modeldata <- namespaceConvertData(model@data, model@name)	
	if (is.null(defaultData)) {
		defaultData <- modeldata
	} 	

	if (!is.null(fitfunction)) {
		if ("data" %in% slotNames(fitfunction)) {
			if(is.na(fitfunction@data) && is.null(modeldata) && !is.null(defaultData)) {
				fitfunction@data <- defaultData@name
			} else if (is.na(fitfunction@data) && !is.null(modeldata)) {
				fitfunction@data <- modeldata@name
			}
		}
		retval <- list(fitfunction)
	} else {
		retval <- list()
	}

	if (length(model@submodels) > 0) {
		submodel_fitfunctions <- lapply(model@submodels, collectFitFunctionsHelper, namespace, defaultData)		
		submodel_fitfunctions <- unlist(submodel_fitfunctions, recursive = FALSE, use.names = FALSE)
		retval <- append(retval, submodel_fitfunctions)
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
