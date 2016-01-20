#
#   Copyright 2007-2016 The OpenMx Project
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

##' imxExtractMethod
##'
##' This is an internal function exported for those people who know
##' what they are doing.
##'
##' @param model model
##' @param index index
imxExtractMethod <- function(model, index) {
	if (is.null(index)) {
		return(NULL)
    }
	if (!(length(index) == 1 && is.character(index))) {
		msg <- paste("The argument to the '$' or '['",
			"operator applied on a MxModel object",
			"must be a single character string")
		stop(msg, call. = FALSE)
	}
	return(namespaceSearch(model, index))
}

##' imxReplaceMethod
##'
##' This is an internal function exported for those people who know
##' what they are doing.
##'
##' @param x the thing
##' @param name name
##' @param value value
imxReplaceMethod <- function(x, name, value) {
	return(namespaceSearchReplace(x, name, value))
}

namespaceSearch <- function(model, name) {
	if (is.na(name) || is.null(name) || identical(name, "")) {
		return(NULL)
	}
	components <- unlist(strsplit(name, imxSeparatorChar, fixed = TRUE))
	if (length(components) == 1) {
		path <- namespaceFindPath(model, name)
		if (is.null(path)) {
			return(namespaceLocalSearch(model, name))
		} else {
			return(namespaceGetModel(model, path))
		}
	} else {
		path <- namespaceFindPath(model, components[[1]])
		if (is.null(path)) {
			return(NULL)
		} else {
			submodel <- namespaceGetModel(model, path)
			return(namespaceLocalSearch(submodel, components[[2]]))			
		}
	}
}

namespaceFindPath <- function(model, modelname) {
	parentname <- model@name
	if (parentname == modelname) {
		return(modelname)
	} else if (length(model@submodels) == 0) {
		return(NULL)
	} else if (modelname %in% names(model@submodels)) {
		return(c(parentname, modelname))
	} else {
		candidates <- lapply(model@submodels, namespaceFindPath, modelname)
		filter <- !(sapply(candidates, is.null))
		ncand <- which(filter)
		ncandlen <- length(ncand)
		if (ncandlen == 0) {
			return(NULL)
		} else if (ncandlen == 1) {
			return(c(parentname, candidates[[ncand]]))
		} else {
			msg <- paste("There are two models with the name",
				omxQuotes(modelname), "found in the model",
				omxQuotes(parentname))
			stop(msg, call. = FALSE)
		}
	}
}

namespaceGetModel <- function(model, path) {
	pathlen <- length(path)
	if (pathlen == 0) {
		stop("An internal error has occured in namespaceSearch")
	} else if (pathlen == 1) {
		return(model)
	} else if (pathlen == 2) {
		return(model@submodels[[path[[2]]]])
	} else {
		remainder <- path[2:pathlen]
		got <- namespaceGetModel(model@submodels[[path[[2]]]], remainder)
		return(got)
	}
}


#
# Check for a named entity within the local model
#
namespaceLocalSearch <- function(model, name) {
	if (name == model@name) {
		return(model)
	}

	if (identical(name, "objective"))  {
		return(list(model@expectation, model@fitfunction))
	}

	if (!is.null(model@data) && name == model@data@name) {
		return(model@data)
	}

	if (!is.null(model@fitfunction) && name == model@fitfunction@name) {
		return(model@fitfunction)
	}

	if (!is.null(model@expectation) && name == model@expectation@name) {
		return(model@expectation)
	}

	result <- model@matrices[[name]]

	if (!is.null(result)) {
		return(result)
	}

	result <- model@algebras[[name]]

	if (!is.null(result)) {
		return(result)
	}

	result <- model@submodels[[name]]

	if (!is.null(result)) {
		return(result)
	}

	result <- model@constraints[[name]]

	if (!is.null(result)) {
		return(result)
	}

	return(NULL)
}

namespaceSearchReplace <- function(model, name, value) {
	if (is.na(name) || identical(name, "")) {
		return(model)
	}
	if (!is.null(value) && !(isS4(value) && ("name" %in% slotNames(value)))) {
		stop(paste("Right-hand side of assignment",
		"operator has illegal value", omxQuotes(value)), call. = FALSE)
	}
	model@.modifiedSinceRun <- TRUE
	components <- unlist(strsplit(name, imxSeparatorChar, fixed = TRUE))
	if (length(components) == 1) {
		path <- namespaceFindPath(model, name)
		if (is.null(path)) {
			return(localNamespaceSearchReplace(model, name, value))
		} else {
			return(namespaceSearchReplaceHelper(model, name, value, path))
		}
	} else {
		path <- namespaceFindPath(model, components[[1]])
		if (is.null(path)) {
			msg <- paste("Could not find the submodel", omxQuotes(components[[1]]),
				"in the interior of model", omxQuotes(model@name))
			stop(msg, call. = FALSE)
		} else {
			return(namespaceSearchReplaceHelper(model, components[[2]], value, path))
		}
	}
}


namespaceSearchReplaceHelper <- function(model, name, value, path) {
	pathlen <- length(path)
	if (pathlen == 0) {
		stop("An internal error has occured in namespaceSearchReplace")
	} else if (pathlen == 1) {
		return(localNamespaceSearchReplace(model, name, value))
	}
	subname <- path[[2]]
	if (pathlen == 2) {
		submodel <- model@submodels[[subname]]
		submodel <- localNamespaceSearchReplace(submodel, name, value)
		model@submodels[[subname]] <- submodel
		return(model)
	} else {
		remainder <- path[2:pathlen]
		model@submodels[[subname]] <- namespaceSearchReplaceHelper(
			model@submodels[[subname]], name, value, remainder)
		return(model)
	}
}

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
	}
	if (identical(name, "objective")) {
		modelname <- model@name
		if (is.null(value)) {
			model <- localNamespaceSearchReplace(model, paste(modelname, "expectation", sep = "."), NULL)
			model <- localNamespaceSearchReplace(model, paste(modelname, "fitfunction", sep = "."), NULL)
			return(model)
		} else if (length(value) == 2) {
			model <- localNamespaceSearchReplace(model, paste(modelname, "expectation", sep = "."), value[[1]])
			model <- localNamespaceSearchReplace(model, paste(modelname, "fitfunction", sep = "."), value[[2]])
			return(model)
		} else {
			msg <- paste("Error in replacing", omxQuotes(name), "in model",
				omxQuotes(model@name), ": value must be either NULL or a (expectation, fit) list")
			stop(msg, call. = FALSE)
		}
	}
	current <- namespaceLocalSearch(model, name)
	if (is.null(current) && is.null(value)) {
		if (!is.na(match(name, c("matrices", "algebras", "constraints",
					 "intervals", "submodels")))) {
			msg <- paste("I'm very sorry, but direct modification of objects inside an mxModel is not",
				     "supported.  The recommended approach to modifying an mxMatrix, mxAlgebra or",
				     "other object within a model is to rebuild the mxModel.  So, for example,",
				     "'model <- mxModel(model, Something)' can be used to replace Something",
				     "inside the model, instead of 'model$Something <- Something' which does not work")
			stop(msg, call. = FALSE)
		}
		return(model)
	}
	if(!is.null(current) && !is.null(value) && 
			!imxSameType(current, value)) {
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
	} else if (is(test,"MxFitFunction")) {
		model@fitfunction <- value
	} else if (is(test,"MxExpectation")) {
		model@expectation <- value
	} else if (is(test,"MxData")) {
		model@data <- value
	} else if (is(test,"MxConstraint")) {
		model@constraints[[name]] <- value
	} else if (is(test, "MxCompute")) {
		model@compute <- value
	} else {
		stop(paste(deparse(test), "is of unknown value for replacement using name",
			name, "in model", model@name), call. = FALSE)
	}
	return(model)
}

