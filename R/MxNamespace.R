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

omxReservedNames <- c('likelihood')

omxVerifyName <- function(name) {
	checkPeriod <- grep('.', name, fixed = TRUE)
	if(length(checkPeriod) > 0) {
		stop(paste("The name", omxQuotes(name),
			"is illegal because it contains a period character."),
			call. = FALSE)
	}
	if (name %in% omxReservedNames) {
		stop(paste("The name", omxQuotes(name),
			"is illegal because it is a reserved name."),
			call. = FALSE)
	}
}

omxIdentifier <- function(namespace, name) {
	return(paste(namespace, name, sep = '.'))
}

omxReverseIdentifier <- function(model, name) {
	components <- unlist(strsplit(name, '.', fixed = TRUE))
	if(length(components) == 1) {
		namespace <- model@name
	} else if (length(components) == 2) {
		namespace <- components[[1]]
		name <- components[[2]]
	} else {
		stop(paste("The reference", omxQuotes(name),
			"has multiple period characters."),
			call. = FALSE)
	}
	return(c(namespace, name))
}

namespaceSearch <- function(model, namespace, name) {
	if (namespace == model@name) {
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
		} else {
			return(fourth)
		}
	} else {
		if (length(model@submodels) == 0) {
			return(NULL)
		}
		recursive <- lapply(model@submodels, 
			function(x) { namespaceSearch(x, namespace, name) })
		target <- which(!sapply(recursive, is.null))
		if (length(target) == 0) {
			return(NULL)
		} else if (length(target) == 1) {
			return(recursive[[target]])
		} else {
			stop(paste("There are two named entities",
				"that matched to the identifier",
				omxQuotes(omxIdentifier(namespace, name))),
				call. = FALSE)
		}
	}
}

namespaceSearchReplace <- function(model, namespace, name, value) {
	if (namespace == model@name) {
		current <- model[[name]]
		if (is.null(current) && is.null(value)) {
			return(model)
		}
		if(name == model@name) {
			stop(paste(omxQuotes(name), 
				"is already used as the name of the model"), call. = FALSE)
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
	} else {
		model@submodels <- lapply(model@submodels,
			function(x) { namespaceSearchReplace(x, namespace, name, value) })
	}
	return(model)
}

