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

flatNamespaceSearch <- function(model, name) {
	if (is.na(name) || is.null(name) || identical(name, "")) {
		return(NULL)
	}
	if (model@name == name) return(model)
	components <- unlist(strsplit(name, imxSeparatorChar, fixed = TRUE))
	if (length(components) == 1) {
		name <- imxIdentifier(model@name, name)
	}
	return(doFlatNamespaceSearch(model, name))
}

flatNamespaceSearchReplace <- function(model, name, value) {
	if (is.na(name) || identical(name, "")) {
		return(model)
	}
	if (!is.null(value) && !(isS4(value) && ("name" %in% slotNames(value)))) {
		stop(paste("Right-hand side of assignment",
		"operator has illegal value", omxQuotes(value)), call. = FALSE)
	}
	if (model@name == name) {
		return(flatNamespaceModelSearchReplace(name, value))
	}
	components <- unlist(strsplit(name, imxSeparatorChar, fixed = TRUE))
	if (length(components) == 1) {
		name <- imxIdentifier(model@name, name)
	}
	return(doFlatNamespaceSearchReplace(model, name, value))
}

#
# Check for a named entity within the flat model
#
doFlatNamespaceSearch <- function(model, name) {

	components <- unlist(strsplit(name, imxSeparatorChar, fixed = TRUE))
	modelname <- components[[1]]
	localname <- components[[2]]

	if (identical(localname, "objective")) {
		expectation <- doFlatNamespaceSearch(model, paste(modelname, "expectation", sep = "."))
		fitfunction <- doFlatNamespaceSearch(model, paste(modelname, "fitfunction", sep = "."))
		return(list(expectation, fitfunction))
	}

	result <- model@matrices[[name]]

	if (!is.null(result)) {
		return(result)
	}

	result <- model@algebras[[name]]

	if (!is.null(result)) {
		return(result)
	} 

	result <- model@constraints[[name]]

	if (!is.null(result)) {
		return(result)
	}

	result <- model@fitfunctions[[name]]

	if (!is.null(result)) {
		return(result)
	}

	result <- model@expectations[[name]]

	if (!is.null(result)) {
		return(result)
	}

	result <- model@datasets[[name]]

	if (!is.null(result)) {
		return(result)
	}

	return(NULL)
}

#
# Provide the list of named entities valid in the namespace
#
flatNamespaceList <- function(model) {
	result <- c(names(model@matrices), 
	names(model@algebras),
	names(model@constraints),
	names(model@fitfunctions),
	names(model@expectations),
	names(model@datasets))
	return(result)
}


#
# First portion of flatNamespaceSearchReplace(model, name, value)
# Replaces the entire flat model if name matches
#
flatNamespaceModelSearchReplace <- function(name, value) {
	if (!(is.null(value) || is(value, "MxFlatModel"))) {
		stop(paste("Replacement for model", omxQuotes(name),
			"is neither NULL nor MxFlatModel object"), call. = FALSE)
	}
	if (!is.null(value)) {
		value@name <- name
	}
	return(value)
}

#
# Replace a named entity within the flat model
#
doFlatNamespaceSearchReplace <- function(model, name, value) {

	components <- unlist(strsplit(name, imxSeparatorChar, fixed = TRUE))
	modelname <- components[[1]]
	localname <- components[[2]]

	if (identical(localname, "objective")) {
		if (is.null(value)) {
			model <- doFlatNamespaceSearchReplace(model, paste(modelname, "expectation", sep = "."), NULL)
			model <- doFlatNamespaceSearchReplace(model, paste(modelname, "fitfunction", sep = "."), NULL)
			return(model)
		} else if (length(value) == 2) {
			model <- doFlatNamespaceSearchReplace(model, paste(modelname, "expectation", sep = "."), value[[1]])
			model <- doFlatNamespaceSearchReplace(model, paste(modelname, "fitfunction", sep = "."), value[[2]])
			return(model)
		} else {
			msg <- paste("Error in replacing", omxQuotes(name), "in model",
				omxQuotes(model@name), ": value must be either NULL or a (expectation, fit) list")
			stop(msg, call. = FALSE)
		}
	}

	current <- doFlatNamespaceSearch(model, name)
	if (is.null(current) && is.null(value)) {
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
	} else if (is(test,"MxFitFunction")) {
		model@fitfunctions[[name]] <- value
	} else if (is(test,"MxExpectation")) {
		model@expectations[[name]] <- value
	} else if (is(test,"MxData")) {
		model@datasets[[name]] <- value
	} else if (is(test,"MxConstraint")) {
		model@constraints[[name]] <- value
	} else {
		stop(paste(test, "is of unknown value for replacement using name",
			name, "in model", model@name), call. = FALSE)
	}
	return(model)
}
