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

namespaceSearch <- function(model, name) {
	if (is.na(name) || is.null(name) || identical(name, "")) {
		return(NULL)
	}
	components <- unlist(strsplit(name, omxSeparatorChar, fixed = TRUE))
	if (length(components) == 1) {
		submodel <- namespaceModelSearch(model, name)
		if (!is.null(submodel)) {
			return(submodel)
		}
	}
	tuple <- omxReverseIdentifier(model, name)
	return(namespaceSearchHelper(model, tuple[[1]], tuple[[2]]))
}

namespaceSearchReplace <- function(model, name, value) {
	if (is.na(name) || identical(name, "")) {
		return(model)
	}
	if (!is.null(value) && !(isS4(value) && ("name" %in% slotNames(value)))) {
		stop(paste("Right-hand side of assignment",
		"operator has illegal value", omxQuotes(value)), call. = FALSE)
	}
	components <- unlist(strsplit(name, omxSeparatorChar, fixed = TRUE))
	if (length(components) == 1) {
		submodel <- namespaceModelSearch(model, name)
		if (!is.null(submodel)) {
			model <- namespaceModelSearchReplace(model, name, value)
			return(model)
		}
	}
	tuple <- omxReverseIdentifier(model, name)
	return(namespaceSearchReplaceHelper(model, tuple[[1]], tuple[[2]], value))
}

#
# First portion of namespaceSearch(model, name)
# Checks for an existing submodel with the same name
#
namespaceModelSearch <- function(model, name) {
	if (model@name == name) return(model)
	results <- lapply(model@submodels, namespaceModelSearch, name)
	names(results) <- NULL
	results <- unlist(results)
	if (length(results) == 0) {
		return(NULL)
	} else if (length(results) == 1) {
		return(results[[1]])
	} else {
		stop(paste("There are two models with the name",
			omxQuotes(name), "found in the model"), call. = FALSE) 
	}
}

#
# Check for a named entity within the local model
#
namespaceLocalSearch <- function(model, name) {
	if (name == model@name) {
		return(model)
	} else if (name %in% names(omxReservedNames) && 
		!is.null(omxReservedNames[[name]]@search)) {
		return(omxReservedNames[[name]]@search(model))
	}
	first <- model@matrices[[name]]
	second <- model@algebras[[name]]
	third <- model@submodels[[name]]
	fourth <- model@constraints[[name]]
	if (!is.null(model@objective) && name == model@objective@name) {
		return(model@objective)
	} else if (!is.null(model@data) && name == model@data@name) {
		return(model@data)
	} else if (!is.null(first)) {
		return(first)
	} else if (!is.null(second)) {
		return(second)
	} else if (!is.null(third)) {
		return(third)
	} else if (!is.null(fourth)) {
		return(fourth)
	}
	return(NULL)
}

#
# The recursive portion of namespaceSearch(model, name)
#
namespaceSearchHelper <- function(model, namespace, name) {
	if (namespace == model@name) {
		return(namespaceLocalSearch(model, name))
	} else {
		if (length(model@submodels) == 0) {
			return(NULL)
		}
		results <- lapply(model@submodels, namespaceSearchHelper,
			namespace, name)
		names(results) <- NULL
		results <- unlist(results)
		if (length(results) == 0) {
			return(NULL)
		} else if (length(results) == 1) {
			return(results[[1]])
		} else {
			stop(paste("There are two named entities",
				"that matched to the identifier",
				omxQuotes(omxIdentifier(namespace, name))),
				call. = FALSE)
		}
	}
}

#
# First portion of namespaceSearchReplace(model, name, value)
# Replaces an existing submodel with the same name
#
namespaceModelSearchReplace <- function(model, name, value) {
	if (!(is.null(value) || is(value, "MxModel"))) {
		stop(paste("Replacement for model", omxQuotes(name),
			"is neither NULL nor MxModel object"), call. = FALSE)
	}
	if (model@name == name) {
		if (!is.null(value)) {
			value@name <- name
		}
		return(value)
	}
	if (length(model@submodels) > 0) {
		submodels <- lapply(model@submodels, namespaceModelSearchReplace, name, value)
		select <- sapply(submodels, is.null)
		model@submodels <- submodels[!select]
	}
	return(model)
}

#
# Replace a named entity within the local model
#
localNamespaceSearchReplace <- function(model, name, value) {
	if (name == model@name) {
		if (!(is.null(value) || is(value, "MxModel"))) {
			stop(paste("Replacement for model", omxQuotes(name),
				"is neither NULL nor MxModel object"), call. = FALSE)
		}
		if (!is.null(value)) {
			value@name <- name
		}
		return(value)
	} else if (name %in% names(omxReservedNames) && 
		!is.null(omxReservedNames[[name]]@replace)) {
		return(omxReservedNames[[name]]@replace(model, value))
	}
	current <- namespaceLocalSearch(model, name)
	if (is.null(current) && is.null(value)) {
		return(model)
	}
	if(!is.null(current) && !is.null(value) && 
			!omxSameType(current, value)) {
		stop(paste("There already exists an object", 
				omxQuotes(name), 
				"in this model of different type"), call. = FALSE)
	}
	if(!is.null(value)) {
		value@name <- name
		test <- value
	} else {
		test <- current
	}
	if (is(test,"MxMatrix")) {
		model@matrices[[name]] <- value
	} else if (is(test,"MxAlgebra")) {
		model@algebras[[name]] <- value
	} else if (is(test,"MxModel")) {
		model@submodels[[name]] <- value	
	} else if (is(test,"MxObjective")) {
		model@objective <- value
	} else if (is(test,"MxData")) {
		model@data <- value
	} else if (is(test,"MxConstraint")) {
		model@constraints[[name]] <- value
	} else {
		stop(paste(test, "is of unknown value for replacement using name",
			name, "in model", model@name), call. = FALSE)
	}
	return(model)
}

#
# The recursive portion of namespaceSearchReplace(model, name, value)
#
namespaceSearchReplaceHelper <- function(model, namespace, name, value) {
	if (namespace == model@name) {
		return(localNamespaceSearchReplace(model, name, value))
	} else {
		submodels <- lapply(model@submodels, 
			namespaceSearchReplaceHelper, namespace, name, value)
		select <- sapply(submodels, is.null)
		model@submodels <- submodels[!select]
	}
	return(model)
}
