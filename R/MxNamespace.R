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

##' imxSeparatorChar
##'
##' The character between the model name and the named entity inside
##' the model.
imxSeparatorChar <- '.'

isNumber <- function(input) {
    match <- grep("^[0-9]+[.]*[0-9]*(L|((E|e)[-+]?[0-9]+))?$", input, perl = TRUE, value = TRUE)
    return(length(match) > 0)
}

explode <- function(string) { 
	strsplit(string, split=character())[[1]] 
}

illegalChars <- "+-!~?:*/^%<>=&|$"
illegalCharsVector <- explode(illegalChars)

availableName <- function(model, namespace, name) {
	return(is.null(model[[name]]) && 
			!(name %in% namespace$parameters) &&
			!(name %in% namespace$values))
}

##' imxVerifyReference
##'
##' This is an internal function exported for those people who know
##' what they are doing.
##'
##' @param reference reference
##' @param stackNumber stackNumber
imxVerifyReference <- function(reference, stackNumber) {
	if (length(reference) != 1) {
		stop(paste("Internal error of call to imxVerifyReference in",
			deparse(width.cutoff = 400L, sys.call(stackNumber - 1))),
		call. = FALSE)
	}
	if (is.na(reference)) {
		return()
	}
    if (isNumber(reference)) {
        stop(paste("The reference", omxQuotes(reference),
            "in", deparse(width.cutoff = 400L, sys.call(stackNumber - 1)), 
			"is illegal because it can be interpreted",
            "as a number"), call. = FALSE)
    }
    if (identical(reference, "")) {
        stop(paste("The reference", omxQuotes(reference),
            "in", deparse(width.cutoff = 400L, sys.call(stackNumber - 1)), 
			"is illegal because references",
            "of zero length are not allowed"), call. = FALSE)
    }
	if (!is.na(reference) && substring(reference, nchar(reference), 
				nchar(reference)) == imxSeparatorChar) {
			stop(paste("The reference", omxQuotes(reference),
				"in", deparse(width.cutoff = 400L, sys.call(stackNumber - 1)), 
				"is illegal because it contains a",
				omxQuotes(imxSeparatorChar), 
				"with either a missing prefix or suffix."),
			call. = FALSE)
	}
	components <- unlist(strsplit(reference, imxSeparatorChar, fixed = TRUE))
	for(i in 1:length(components)) {
		component <- components[[i]]
		if (nchar(component) == 0) {
			stop(paste("The reference", omxQuotes(reference),
				"in", deparse(width.cutoff = 400L, sys.call(stackNumber - 1)), 
				"is illegal because it contains a",
				omxQuotes(imxSeparatorChar), 
				"with either a missing prefix or suffix."),
			call. = FALSE)
		}
	}
	if (length(components) == 2) {
		if ((components[[1]] != "data") && !hasSquareBrackets(reference)) {
			stop(paste("The reference", omxQuotes(reference),
				"is illegal because it contains the",
				omxQuotes(imxSeparatorChar), "character in", 
				deparse(width.cutoff = 400L, sys.call(stackNumber - 1)),
				". To write a definition variable use", 
				omxQuotes(paste("data", components[[2]], sep = "."))),
			call. = FALSE)
		}
	} else if (length(components) == 3) {
		if ((components[[2]] != "data")) {
			stop(paste("The reference", omxQuotes(reference),
				"is illegal because it contains the",
				omxQuotes(imxSeparatorChar), "character",
				"but it is not a valid definition variable in", 
				deparse(width.cutoff = 400L, sys.call(stackNumber - 1))),
			call. = FALSE)
		}
	} else if (length(components) > 3) {
			stop(paste("The reference", omxQuotes(reference),
				"is illegal because it contains the",
				omxQuotes(imxSeparatorChar), "character",
				"but it is not a valid definition variable in", 
				deparse(width.cutoff = 400L, sys.call(stackNumber - 1))),
			call. = FALSE)
	}
	badCharacterMatch <- illegalCharsVector %in% explode(reference)
	if(any(badCharacterMatch)) {
		badChars <- illegalCharsVector[badCharacterMatch]
		if(length(badChars) == 1) {
			stop(paste("The reference", omxQuotes(reference),
				"in", deparse(width.cutoff = 400L, sys.call(stackNumber - 1)), 
				"is illegal because it contains the",
				omxQuotes(badChars[[1]]), 
				"character."),
			call. = FALSE)		
		} else {
			stop(paste("The reference", omxQuotes(reference),
				"in", deparse(width.cutoff = 400L, sys.call(stackNumber - 1)), 
				"is illegal because it contains the characters",
				paste(omxQuotes(badChars), '.', sep = '')),
			call. = FALSE)		
		}
	}
	leftSquareBracket <- length(grep("[", reference, fixed=TRUE)) > 0
	rightSquareBracket <- length(grep("]", reference, fixed=TRUE)) > 0
	comma <- length(grep(",", reference, fixed=TRUE)) > 0
	if ((leftSquareBracket && rightSquareBracket && comma) ||
		(!leftSquareBracket && !rightSquareBracket && !comma)) {
	} else {
			stop(paste("The reference", omxQuotes(reference),
				"in", deparse(width.cutoff = 400L, sys.call(stackNumber - 1)), 
				"is illegal because it is",
				"a partial square-bracket reference."),
			call. = FALSE)
	}
}

#' mxMakeNames
#'
#' Adjust a character vector so that it can be used as MxMatrix column
#' or row names. OpenMx is (much) more restrictive than base R's make.names.
#'
#' @param names a character vector
#' @param unique whether the pass the result through \link[base]{make.unique}
#' @seealso
#' \link[base]{make.names}
#' @examples
#' demo <- c("", "103", "data", "foo.bar[3,2]", "+!", "!+")
#' mxMakeNames(demo, unique=TRUE)
mxMakeNames <- function(names, unique = FALSE) {
	names <- gsub("\\s", "", names, perl=TRUE)
	names[nchar(names) == 0] <- 'i'
	names <- sapply(names, function(str) {
		if (isNumber(str)) {
			str <- paste("X",str,sep="")
		}
		str
	})
	names[!is.na(match(names, imxReservedNames))] <- "reserved"
	broadlyIllegal <- paste("[\\Q", illegalChars, ".[]\\E]", sep="")
	names <- gsub(broadlyIllegal, "x", names, perl=TRUE)

	if (unique) {
		names <- make.unique(names, sep="")
	}
	names
}

##' imxVerifyName
##'
##' This is an internal function exported for those people who know
##' what they are doing.
##'
##' @param name name
##' @param stackNumber stackNumber
imxVerifyName <- function(name, stackNumber) {
    if (length(name) == 0) {
        stop(paste("The empty character vector is an invalid name in", 
			deparse(width.cutoff = 400L, sys.call(stackNumber - 1))),
            call. = FALSE)
    }
    if (length(name) > 1) {
        stop(paste("The 'name' argument must be a single character argument in", 
			deparse(width.cutoff = 400L, sys.call(stackNumber - 1))),
            call. = FALSE)
    }
    if (identical(name, "")) {
        stop(paste("The empty string is an invalid name in", 
			deparse(width.cutoff = 400L, sys.call(stackNumber - 1))), call. = FALSE)
    }
    if (isNumber(name)) {
        stop(paste("The name", omxQuotes(name),
            "is illegal because it can be interpreted",
            "as a number in", 
			deparse(width.cutoff = 400L, sys.call(stackNumber - 1))), call. = FALSE)
    }
	components <- unlist(strsplit(name, imxSeparatorChar, fixed = TRUE))
	if (length(components) == 2) {
		if (components[[1]] != "data") {
			stop(paste("The name", omxQuotes(name),
				"is illegal because it contains the",
				omxQuotes(imxSeparatorChar), "character in", 
				deparse(width.cutoff = 400L, sys.call(stackNumber - 1))),
			call. = FALSE)			
		}
	} else if (length(components) > 2) {
			stop(paste("The name", omxQuotes(name),
				"is illegal because it contains multiple",
				omxQuotes(imxSeparatorChar), "characters in", 
				deparse(width.cutoff = 400L, sys.call(stackNumber - 1))),
			call. = FALSE)
	}
	if (name %in% imxReservedNames) {
		stop(paste("The name", omxQuotes(name),
			"is illegal because it is a reserved name in", 
			deparse(width.cutoff = 400L, sys.call(stackNumber - 1))),
			call. = FALSE)
	}
	components <- unlist(strsplit(name, '[', fixed = TRUE))
	if (length(components) > 1) {
		stop(paste("The name", omxQuotes(name),
			"is illegal because it contains the '[' character in", 
			deparse(width.cutoff = 400L, sys.call(stackNumber - 1))),
			call. = FALSE)
	}
	components <- unlist(strsplit(name, ']', fixed = TRUE))
	if (length(components) > 1) {
		stop(paste("The name", omxQuotes(name),
			"is illegal because it contains the ']' character in", 
			deparse(width.cutoff = 400L, sys.call(stackNumber - 1))),
			call. = FALSE)
	}
	badCharacterMatch <- illegalCharsVector %in% explode(name)
	if(any(badCharacterMatch)) {
		badChars <- illegalCharsVector[badCharacterMatch]
		if(length(badChars) == 1) {
			stop(paste("The name", omxQuotes(name),
				"is illegal because it contains the",
				omxQuotes(badChars[[1]]), 
				"character in", 
				deparse(width.cutoff = 400L, sys.call(stackNumber - 1))),
			call. = FALSE)		
		} else {
			stop(paste("The name", omxQuotes(name),
				"is illegal because it contains the characters",
				paste(omxQuotes(badChars), sep = ''), "in", 
				deparse(width.cutoff = 400L, sys.call(stackNumber - 1))),
			call. = FALSE)
		}
	}	
	
}

##' imxIsDefinitionVariable
##'
##' This is an internal function exported for those people who know
##' what they are doing.
##'
##' @param name name
imxIsDefinitionVariable <- function(name) {
	if (is.na(name)) {
		return(FALSE)
	}
	components <- unlist(strsplit(name, imxSeparatorChar, fixed = TRUE))
	if (length(components) == 2 && components[[1]] == 'data') {
		return(TRUE)
	} else if (length(components) > 2 && components[[2]] == 'data') {
		return(TRUE)
	} else {
		return(FALSE)
	}
}

isLocalDefinitionVariable <- function(name) {
	components <- unlist(strsplit(name, imxSeparatorChar, fixed = TRUE))
	if (length(components) == 2 && components[[1]] == 'data') {
		return(TRUE)
	} else {
		return(FALSE)
	}
}


##' imxHasDefinitionVariable
##'
##' This is an internal function exported for those people who know
##' what they are doing.  This function checks if a model (or its
##' submodels) has at least one definition variable.
##'
##' @param model model
imxHasDefinitionVariable <- function(model) {
	# Check submodels for defvar
	if(length(model$submodels) > 0){
		for(i in 1:length(model@submodels)){
			attempt <- sapply(model@submodels[[i]], imxHasDefinitionVariable)
			if(any(attempt)){
				return(TRUE)
			}
		}
	}
	
	# Check if the model has data
	if(length(model@data) == 0){
		return(FALSE)
	}
	
	# Check the matrices for defvar
	if(length(model@matrices) > 0){
		for(i in 1:length(model@matrices)){
			attempt <- sapply(model@matrices[[i]]$labels, imxIsDefinitionVariable)
			if(any(attempt)){
				return(TRUE)
			}
		}
	}
	
	# Check the algebras for defvar
	if(length(model@algebras) > 0){
		for(i in 1:length(model@algebras)){
			attempt <- sapply(as.character(model@algebras[[i]]$formula), imxIsDefinitionVariable)
			if(any(attempt)){
				return(TRUE)
			}
		}
	}
	
	# All checks find nothing, return FALSE
	return(FALSE)
}

##' imxIdentifier
##'
##' This is an internal function exported for those people who know
##' what they are doing.
##'
##' @param namespace namespace
##' @param name name
imxIdentifier <- function(namespace, name) {
	if (is.null(name)) return(NULL)
	return(paste(namespace, name, sep = imxSeparatorChar))
}

simplifyName <- function(flatName, modelName) {
	components <- unlist(strsplit(flatName, imxSeparatorChar, fixed = TRUE))
	if (length(components) == 2 && components[[1]] == modelName) {
		return(components[[2]])
	} else if (length(components) == 3 && components[[1]] == modelName) {
		return(paste(components[[2]], components[[3]], sep = '.'))
	} else {
		return(flatName)
	}
}

##' imxReverseIdentifier
##'
##' This is an internal function exported for those people who know
##' what they are doing.
##'
##' @param model model
##' @param name name
imxReverseIdentifier <- function(model, name) {
	components <- unlist(strsplit(name, imxSeparatorChar, fixed = TRUE))
	if(length(components) < 2) {
		namespace <- model@name
	} else if (length(components) == 2) {
		if (components[[1]] == 'data') {
			namespace <- model@name
		} else {
			namespace <- components[[1]]
			name <- components[[2]]
		}
	} else if (length(components) == 3) {
		if (components[[2]] == 'data') {
			namespace <- components[[1]]
			name <- paste(components[[2]], components[[3]], sep = '.')
		} else {
			stop(paste("The reference", omxQuotes(name),
				"has multiple separator characters."),
				call. = FALSE)
		}
	} else {
		stop(paste("The reference", omxQuotes(name),
			"has multiple separator characters."),
			call. = FALSE)
	}
	return(c(namespace, name))
}

##' imxGenerateNamespace
##' 
##' This is an internal function exported for those people who know
##' what they are doing.
##'
##' @param model model
imxGenerateNamespace <- function(model) {
	entities <- list()
	result <- generateLocalNamespace(model)
	entities[[model@name]] <- result[[1]]
	parameters <- result[[2]]
	values <- result[[3]]
	results <- omxLapply(model@submodels, imxGenerateNamespace)
	if (length(results) > 0) {
		for (i in 1:length(results)) {
			subentities <- results[[i]][['entities']]
			subnames <- names(subentities)
			for (j in 1:length(subentities)) {
				if (!is.null(entities[[subnames[[j]]]])) {
					stop(namespaceErrorMessage(subnames[[j]]), call. = FALSE)
				}
				entities[[subnames[[j]]]] <- subentities[[j]]
			}
			parameters <- union(parameters, results[[i]][['parameters']])
            values <- union(values, results[[i]][['values']])
		}
	}
	return(list('entities' = entities, 'parameters' = parameters, 'values' = values))
}

getEntities <- function(namespace) {
	return(namespace[['entities']])	
}

getParameters <- function(namespace) {
	return(namespace[['parameters']])	
}

getValues <- function(namespace) {
	return(namespace[['values']])	
}

generateLocalNamespace <- function(model) {
	thisEntities <- character()
	thisEntities <- namespaceGetEntities(model, "matrices", thisEntities)
	thisEntities <- namespaceGetEntities(model, "algebras", thisEntities)
	thisEntities <- namespaceGetEntities(model, "submodels", thisEntities)
	thisEntities <- namespaceGetEntities(model, "constraints", thisEntities)
	fitfunction <- model@fitfunction
	expectation <- model@expectation
	data <- model@data

	if (!is.null(fitfunction) && (fitfunction@name %in% thisEntities)) {
		stop(namespaceErrorMessage(fitfunction@name), call. = FALSE)
	} else if (!is.null(fitfunction)) {
		thisEntities <- c(thisEntities, fitfunction@name)
		thisEntities <- c(thisEntities, genericFitNewEntities(fitfunction))
	}

	if (!is.null(expectation) && (expectation@name %in% thisEntities)) {
		stop(namespaceErrorMessage(expectation@name), call. = FALSE)
	} else if (!is.null(expectation)) {
		thisEntities <- c(thisEntities, expectation@name)
	}

	if (!is.null(data) && (data@name %in% thisEntities)) {
		stop(namespaceErrorMessage(data@name), call. = FALSE)
	} else if (!is.null(data)) {
		thisEntities <- c(thisEntities, data@name)
	}

	thisParameters <- namespaceGetParameters(model)
    thisValues <- namespaceGetValues(model)
	return(list('entities' = thisEntities, 'parameters' = thisParameters, 'values' = thisValues))
}

namespaceGetParameters <- function(model) {
	parameters <- sapply(model@matrices, function(x) {
			labels <- x@labels
			labels <- unique(labels[!is.na(labels) & x@free])
			return(labels)
		})
	parameters <- unlist(parameters)
	names(parameters) <- NULL
	return(parameters)
}

namespaceGetValues <- function(model) {
	values <- sapply(model@matrices, function(x) {
			labels <- x@labels
			labels <- unique(labels[!is.na(labels) & !x@free])
			defVars <- sapply(labels, imxIsDefinitionVariable)
			labels <- labels[!defVars]
			subs <- sapply(labels, hasSquareBrackets)
			labels <- labels[!subs]
			return(labels)
		})
	values <- unlist(values)
	names(values) <- NULL
	return(values)
}

namespaceGetEntities <- function(model, slotname, thisEntities) {
	entities <- slot(model, slotname)
	entityNames <- names(entities)
	checkNameAlignment(entityNames, imxExtractNames(entities))
	entityIntersect <- intersect(entityNames, thisEntities)
	if (length(entityIntersect) > 0) {
		stop(namespaceErrorMessage(entityIntersect), call. = FALSE)
	}
	thisEntities <- c(thisEntities, entityNames)
	return(thisEntities)
}

namespaceErrorMessage <- function(rlist) {
	if (length(rlist) == 1) {
		return(paste("The name", omxQuotes(rlist), 
		"appears more than once in the model\n"))
	} else {
		return(paste("The names", omxQuotes(rlist), 
		"appear more than once in the model\n"))
	}
}

##' imxExtractNames
##'
##' This is an internal function exported for those people who know
##' what they are doing.
##'
##' @param lst lst
imxExtractNames <- function(lst) {
	if (length(lst) == 0) {
		return(character())
	} else {
		return(sapply(lst, slot, 'name'))
	}
}

##' imxExtractReferences
##'
##' This is an internal function exported for those people who know
##' what they are doing.
##'
##' @param lst lst
imxExtractReferences <- function(lst) {
	if (length(lst) == 0) {
		return(character())
	} else {
		return(sapply(lst, slot, 'reference'))
	}
}

checkNameAlignment <- function(lst1, lst2) {
	if(any(lst1 != lst2)) {
		stop(paste("The names for entities", omxQuotes(lst1[lst1 != lst2]),
			"do not match their designations in the model"), call. = FALSE)
	}
}

##' omxCheckNamespace
##'
##' This is an internal function exported for those people who know
##' what they are doing.
##'
##' @param model model
##' @param namespace namespace
##' @details This function checks that the named entities in the model are valid.
omxCheckNamespace <- function(model, namespace) {
	checkNamespaceHelper(model, model, namespace)
}

checkNamespaceHelper <- function(model, topmodel, namespace) {
	lapply(model@matrices, checkNamespaceMatrix, model, namespace)
	lapply(model@algebras, checkNamespaceAlgebra, model, namespace)
	lapply(model@constraints, checkNamespaceConstraint, model, namespace)
	lapply(imxDependentModels(model), checkNamespaceHelper, topmodel, namespace)
	allEntities <- unique(unlist(namespace$entities))
	overlap <- intersect(allEntities, namespace$parameters)
	if (length(overlap) > 0) {
		item1 = overlap[1]
		stop(paste("In model", omxQuotes(model@name),
			"the following are both named",
			"entities and free parameters:",
			omxQuotes(overlap),
			"\nIf you are trying to set a path using an mxAlgebra, then",
			"refer to the Algebra with square-bracket notation.",
			"\ni.,e, instead of labels=\"", omxQuotes(overlap), "\"",
			"use: labels=\"", omxQuotes(overlap), "[1,1]\""), call. = FALSE)
	}
	overlap <- intersect(allEntities, namespace$values)
	if (length(overlap) > 0) {
		stop(paste("In model", omxQuotes(model@name),
			"the following are both named",
			"entities and fixed parameters:",
			omxQuotes(overlap)), call. = FALSE)
	}
	overlap <- intersect(namespace$parameters, namespace$values)
	if (length(overlap) > 0) {
		select <- overlap[[1]]
		freelocation <- imxLocateLabel(select, topmodel, TRUE)
		fixedlocation <- imxLocateLabel(select, topmodel, FALSE)
		stop(paste("In model", omxQuotes(topmodel@name),
			"the name", omxQuotes(select),
			"is used as a free parameter in", omxQuotes(freelocation),
			"and as a fixed parameter in", omxQuotes(fixedlocation)), call. = FALSE)
	}
}

##' imxLocateLabel
##'
##' This is an internal function exported for those people who know
##' what they are doing.
##'
##' @param label label
##' @param model model
##' @param parameter parameter
imxLocateLabel <- function(label, model, parameter) {
	if ((length(label) != 1) || !is.character(label) || is.na(label)) {
		stop("'label' must be a character value")
	}
	if ((length(model) != 1) || !is(model, "MxModel")) {
		stop("'model' must be a MxModel object")
	}
	if ((length(parameter) != 1) || !is.logical(parameter) || is.na(parameter)) {
		stop("'parameter' must be either TRUE or FALSE")
	}
	values <- sapply(model@matrices, locateLabelHelper, model, label, parameter)	
	values <- setdiff(values, '')
	children <- lapply(model@submodels, imxLocateLabel, label = label, parameter = parameter)
	children <- unlist(children, recursive = TRUE)
	retval <- as.character(c(values, children))
	return(retval)
}

locateLabelHelper <- function(matrix, model, label, parameter) {
	result <- (matrix@free == parameter) & !is.na(matrix@labels) & (matrix@labels == label)
	if(any(result)) return(imxIdentifier(model@name, matrix@name))
	else return('')
}

legalGlobalReference <- function(name) {
	if(exists(name, envir = globalenv())) {
		value <- get(name, envir = globalenv())
		return(is.numeric(value) || is.character(value))
	}
	return(FALSE);
}

checkNamespaceIdentifier <- function(identifier, model, entity, namespace) {
	entities <- getEntities(namespace)
	parameters <- getParameters(namespace)
    values <- getValues(namespace)
	identifier <- imxReverseIdentifier(model, identifier)
	space <- identifier[[1]]
	name <- identifier[[2]]
	if ( !(name %in% entities[[space]]) &&
		 !(imxIsDefinitionVariable(name)) &&
		 !(hasSquareBrackets(name)) &&
		 !(name %in% parameters) &&
         !(name %in% values) &&
         !(legalGlobalReference(name)) &&
		 !(name %in% imxReservedNames)) {
		stop(paste("Unknown reference", 
			omxQuotes(simplifyName(imxIdentifier(space, name), model@name)),
			"detected in the entity", omxQuotes(entity),
			"in model", omxQuotes(model@name)), call. = FALSE)
	}
}

checkNamespaceAlgebra <- function(algebra, model, namespace) {
	formula <- algebra@formula
	checkNamespaceFormula(formula, model, algebra, namespace)
}

checkNamespaceFormula <- function(formula, model, entity, namespace) {
	if (length(formula) == 1) {
        if (is.numeric(formula)) {
        } else if (identical(as.character(formula), "")) {
        } else {
    		checkNamespaceIdentifier(as.character(formula), model, entity@name, namespace)
        }
	} else if (length(formula) == 4 && identical(as.character(formula[1]), '[')) {
			checkNamespaceFormula(formula[[2]], model, entity, namespace)
			checkNamespaceFormulaAllowCharacters(formula[[3]], model, entity, namespace)
			checkNamespaceFormulaAllowCharacters(formula[[4]], model, entity, namespace)
	} else {
		for (i in 2:length(formula)) {
			checkNamespaceFormula(formula[[i]], model, entity, namespace)
		}
	}
}

checkNamespaceFormulaAllowCharacters <- function(formula, model, entity, namespace) {
	if (length(formula) == 1) {
        if (is.numeric(formula)) {
        } else if (is.character(formula)) {
        } else if (identical(as.character(formula), "")) {
        } else {
    		checkNamespaceIdentifier(as.character(formula), model, entity@name, namespace)
        }
	} else if (length(formula) == 4 && identical(as.character(formula[1]), '[')) {
			checkNamespaceFormula(formula[[2]], model, entity, namespace)
			checkNamespaceFormulaAllowCharacters(formula[[3]], model, entity, namespace)
			checkNamespaceFormulaAllowCharacters(formula[[4]], model, entity, namespace)
	} else {
		for (i in 2:length(formula)) {
			checkNamespaceFormulaAllowCharacters(formula[[i]], model, entity, namespace)
		}
	}
}


checkNamespaceConstraint <- function(constraint, model, namespace) {
	formula <- constraint@formula
	checkNamespaceFormula(formula, model, constraint, namespace)
}

checkNamespaceMatrix <- function(matrix, model, namespace) {
	labels <- matrix@labels
	notNAlabels <- labels[!is.na(labels) & matrix@free]
	lapply(notNAlabels, function(x) { checkNamespaceIdentifier(x, model, matrix@name, namespace) })
}

##' imxConvertSubstitution
##'
##' This is an internal function exported for those people who know
##' what they are doing.
##'
##' @param substitution substitution
##' @param modelname modelname
##' @param namespace namespace
imxConvertSubstitution <- function(substitution, modelname, namespace) {
	pieces <- splitSubstitution(substitution)
	identifier <- imxConvertIdentifier(pieces[[1]], modelname, namespace)
	result <- paste(identifier, '[', pieces[[2]], ',', pieces[[3]], ']', sep = '')
	return(result)
}

##' imxConvertIdentifier
##' 
##' This is an internal function exported for those people who know
##' what they are doing.
##'
##' @param identifiers identifiers
##' @param modelname modelname
##' @param namespace namespace
##' @param strict strict
imxConvertIdentifier <- function(identifiers, modelname, namespace, strict=FALSE) {
	if (length(identifiers) == 0) return(identifiers)
	identifiers <- as.character(identifiers)
	if (all(is.na(identifiers))) return(identifiers)
	vapply(identifiers, function(identifier) {
		isLocalEntity <- identifier %in% namespace$entities[[modelname]]
		if (isLocalEntity) {
			return(imxIdentifier(modelname, identifier))
		} else if (isLocalDefinitionVariable(identifier)) {
			return(imxIdentifier(modelname, identifier))
		} else {
			if (strict && length(unlist(strsplit(identifier, imxSeparatorChar, fixed = TRUE))) == 1) {
				stop(paste("Identifier", omxQuotes(identifier), "refers to what?"))
			}
			return(identifier)
		}
	}, "", USE.NAMES=FALSE)
}

getModelName <- function(object) {
	return(unlist(strsplit(object@name, imxSeparatorChar, fixed = TRUE))[[1]])
}

getModelNameString <- function(name) {
	return(unlist(strsplit(name, imxSeparatorChar, fixed = TRUE))[[1]])
}

##' imxConvertLabel
##'
##' This is an internal function exported for those people who know
##' what they are doing.
##'
##' @param label label
##' @param modelname modelname
##' @param dataname dataname
##' @param namespace namespace
imxConvertLabel <- function(label, modelname, dataname, namespace) {
	if (hasSquareBrackets(label)) {
		components <- splitSubstitution(label)
		identifier <- imxConvertLabel(components[[1]], modelname, dataname, namespace)
		row <- imxConvertLabel(components[[2]], modelname, dataname, namespace)
		col <- imxConvertLabel(components[[3]], modelname, dataname, namespace)
		results <- paste(identifier, '[', row, ',', col, ']', sep = '')
		return(results)
	}
	components <- unlist(strsplit(label, imxSeparatorChar, fixed = TRUE))
	if (imxIsDefinitionVariable(label)) {
		if (length(components) == 3) {
			return(label)
		} else if (is.null(dataname)) {
			stop(paste("A definition variable", omxQuotes(label),
				"has been declared in model", omxQuotes(modelname),
				"that does not contain a data set"), call. = FALSE)
		} else {
			datasource <- unlist(strsplit(dataname, imxSeparatorChar, fixed = TRUE))[[1]]
			return(imxIdentifier(datasource, label))
		}

	}
	return(imxConvertIdentifier(label, modelname, namespace))
}

qualifyNamesMatrix <- function(matrix, modelname, dataname, namespace) {
	matrix@name <- imxIdentifier(modelname, matrix@name)
	free <- matrix@free
	labels <- matrix@labels
	select <- (!free) & (!is.na(labels))
	if (any(select)) {
		refNames <- labels[select]
		matrix@labels[select] <- sapply(refNames, imxConvertLabel, modelname, dataname, namespace)
	}
	return(matrix)
}

qualifyNamesAlgebra <- function(algebra, modelname, namespace) {
	algebra@name <- imxIdentifier(modelname, algebra@name)
	algebra@formula <- qualifyNamesFormula(algebra@formula, modelname, namespace)
	return(algebra)
}

qualifyNamesFormula <- function(formula, modelname, namespace) {
	if (length(formula) == 1) {
        if (is.symbol(formula) && 
            (as.character(formula) %in% namespace$parameters || 
             as.character(formula) %in% namespace$values)) {
        } else if (is.numeric(formula)) {
        } else if (identical(as.character(formula), "")) {
        } else {
            result <- imxConvertIdentifier(formula, modelname, namespace)
            if (is.symbol(formula) && is.character(result)) {
                formula <- as.symbol(result)
            } else {
                formula <- result
            }
        }
	} else {
		for (i in 2:length(formula)) {
			formula[[i]] <- qualifyNamesFormula(formula[[i]], modelname, namespace)
		}
	}
	return(formula)
}

qualifyNamesConstraint <- function(constraint, modelname, namespace) {
	constraint@name <- imxIdentifier(modelname, constraint@name)
	constraint@formula <- qualifyNamesFormula(constraint@formula, modelname, namespace)
	return(constraint)
}

qualifyNamesInterval <- function(interval, modelname, namespace) {
	interval@reference <- imxConvertLabel(interval@reference, modelname, NULL, namespace)
	return(interval)
}

safeQualifyNames <- function(obj, modelname, namespace) {
	if (!is.null(obj)) {
		obj <- qualifyNames(obj, modelname, namespace)
	}
	obj
}

qualifyNamesData <- function(data, modelname) {
	if (!is.null(data)) {
		data@name <- imxIdentifier(modelname, data@name)
	}
	return(data)
}
