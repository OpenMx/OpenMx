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
	if (name %in% names(imxReservedNames) && 
		!is.null(imxReservedNames[[name]]@search)) {
		return(imxReservedNames[[name]]@search(model))
	}
	first <- model@matrices[[name]]
	second <- model@algebras[[name]]
	third <- model@constraints[[name]]
	fourth <- model@objectives[[name]]
	fifth <- model@datasets[[name]]
	if (!is.null(first)) {
		return(first)
	} else if (!is.null(second)) {
		return(second)
	} else if (!is.null(third)) {
		return(third)
	} else if (!is.null(fourth)) {
		return(fourth)
	} else if (!is.null(fifth)) {
		return(fifth)
	}
	return(NULL)
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
	if (name %in% names(imxReservedNames) && 
		!is.null(imxReservedNames[[name]]@replace)) {
		return(imxReservedNames[[name]]@replace(model, value))
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
	} else if (is(test,"MxObjective")) {
		model@objectives[[name]] <- value
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
