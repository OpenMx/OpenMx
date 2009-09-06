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

omxSeparatorChar <- '.'

isNumber <- function(input) {
    match <- grep("^[0-9]+[.]*[0-9]*(L|((E|e)[-+]?[0-9]+))?$", input, perl = TRUE, value = TRUE)
    return(length(match) > 0)
}

omxVerifyReference <- function(reference) {
    if (isNumber(reference)) {
        stop(paste("The reference", omxQuotes(reference),
            "is illegal because it can be interpreted",
            "as a number"), call. = FALSE)
    }
	if (!is.na(reference) && substring(reference, nchar(reference), 
				nchar(reference)) == omxSeparatorChar) {
			stop(paste("The reference", omxQuotes(reference),
				"is illegal because it contains a",
				omxQuotes(omxSeparatorChar), 
				"with either a missing prefix or suffix."),
			call. = FALSE)
	}
	components <- unlist(strsplit(reference, omxSeparatorChar, fixed = TRUE))
	for(i in 1:length(components)) {
		component <- components[[i]]
		if (nchar(component) == 0) {
			stop(paste("The reference", omxQuotes(reference),
				"is illegal because it contains a",
				omxQuotes(omxSeparatorChar), 
				"with either a missing prefix or suffix."),
			call. = FALSE)
		}
	}
}

omxVerifyName <- function(name) {
    if (isNumber(name)) {
        stop(paste("The name", omxQuotes(name),
            "is illegal because it can be interpreted",
            "as a number"), call. = FALSE)
    }
	components <- unlist(strsplit(name, omxSeparatorChar, fixed = TRUE))
	if (length(components) == 2) {
		if (components[[1]] != "data") {
			stop(paste("The name", omxQuotes(name),
				"is illegal because it contains the",
				omxQuotes(omxSeparatorChar), "character"),
			call. = FALSE)			
		}
	} else if (length(components) > 2) {
			stop(paste("The name", omxQuotes(name),
				"is illegal because it contains multiple",
				omxQuotes(omxSeparatorChar), "characters"),
			call. = FALSE)
	}
	if (name %in% names(omxReservedNames)) {
		stop(paste("The name", omxQuotes(name),
			"is illegal because it is a reserved name"),
			call. = FALSE)
	}
}

omxIsDefinitionVariable <- function(name) {
	components <- unlist(strsplit(name, omxSeparatorChar, fixed = TRUE))
	if (length(components) == 2 && components[[1]] == 'data') {
		return(TRUE)
	} else if (length(components) > 2 && components[[2]] == 'data') {
		return(TRUE)
	} else {
		return(FALSE)
	}
}

omxIdentifier <- function(namespace, name) {
	return(paste(namespace, name, sep = omxSeparatorChar))
}

simplifyName <- function(flatName, modelName) {
	components <- unlist(strsplit(flatName, omxSeparatorChar, fixed = TRUE))
	if (length(components) == 2 && components[[1]] == modelName) {
		return(components[[2]])
	} else {
		return(flatName)
	}
}

omxReverseIdentifier <- function(model, name) {
	components <- unlist(strsplit(name, omxSeparatorChar, fixed = TRUE))
	if(length(components) == 1) {
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

omxGenerateNamespace <- function(model) {
	entities <- list()
	result <- generateLocalNamespace(model)
	entities[[model@name]] <- result[[1]]
	parameters <- result[[2]]
    values <- result[[3]]
	results <- lapply(model@submodels, omxGenerateNamespace)
	if (length(results) > 0) {
		for (i in 1:length(results)) {
			submodel <- model@submodels[[i]]
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
	objective <- model@objective
	data <- model@data
	if (!is.null(objective) && (objective@name %in% thisEntities)) {
		stop(namespaceErrorMessage(objective@name), call. = FALSE)
	} else if (!is.null(objective)) {
		thisEntities <- c(thisEntities, objective@name)
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
			return(labels)
		})
	values <- unlist(values)
	names(values) <- NULL
	return(values)
}

namespaceGetEntities <- function(model, slotname, thisEntities) {
	entities <- slot(model, slotname)
	entityNames <- names(entities)
	checkNameAlignment(entityNames, entityExtractNames(entities))
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

entityExtractNames <- function(lst) {
	lapply(lst, function(x) { x@name } )	
}

checkNameAlignment <- function(lst1, lst2) {
	if(any(lst1 != lst2)) {
		stop(paste("The names for entities", omxQuotes(lst1[lst1 != lst2]),
			"do not match their designations in the model"), call. = FALSE)
	}
}

omxCheckNamespace <- function(model, namespace) {
	lapply(model@matrices, function(x) { checkNamespaceMatrix(x, model, namespace) })
	lapply(model@algebras, function(x) { checkNamespaceAlgebra(x, model, namespace) })
	lapply(model@constraints, function(x) { checkNamespaceConstraint(x, model, namespace) })
	lapply(model@submodels, function(x) { omxCheckNamespace(x, namespace) })
	allEntities <- unique(unlist(namespace$entities))
	overlap <- intersect(allEntities, namespace$parameters)
	if (length(overlap) > 0) {
		stop(paste("In model", omxQuotes(model@name),
			"the following are both named",
			"entities and free parameters:",
			omxQuotes(overlap)), call. = FALSE)
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
		stop(paste("In model", omxQuotes(model@name),
			"the following are both free",
			"and fixed parameters:",
			omxQuotes(overlap)), call. = FALSE)
	}
}

checkNamespaceIdentifier <- function(identifier, model, namespace) {
	entities <- getEntities(namespace)
	parameters <- getParameters(namespace)
    values <- getValues(namespace)
	identifier <- omxReverseIdentifier(model, identifier)
	space <- identifier[[1]]
	name <- identifier[[2]]
	if ( !(name %in% entities[[space]]) &&
		 !(omxIsDefinitionVariable(name)) &&
		 !(name %in% parameters) &&
         !(name %in% values) &&
		 !(name %in% names(omxReservedNames))) {
		stop(paste("Unknown reference: ", 
			omxQuotes(omxIdentifier(space, name))), call. = FALSE)
	}
}

checkNamespaceAlgebra <- function(algebra, model, namespace) {
	formula <- algebra@formula
	checkNamespaceFormula(formula, model, namespace)
}

checkNamespaceFormula <- function(formula, model, namespace) {
	if (length(formula) == 1) {
        if (is.numeric(formula)) {
        } else {
    		checkNamespaceIdentifier(as.character(formula), model, namespace)
        }
	} else {
		for (i in 2:length(formula)) {
			checkNamespaceFormula(formula[[i]], model, namespace)
		}
	}
}

checkNamespaceConstraint <- function(constraint, model, namespace) {
	checkNamespaceIdentifier(constraint@alg1, model, namespace)
	checkNamespaceIdentifier(constraint@alg2, model, namespace)
}

checkNamespaceMatrix <- function(matrix, model, namespace) {
	labels <- matrix@labels
	notNAlabels <- labels[!is.na(labels) & matrix@free]
	lapply(notNAlabels, function(x) { checkNamespaceIdentifier(x, model, namespace) })
}

omxConvertIdentifier <- function(identifier, modelname, namespace) {
    isLocalEntity <- as.character(identifier) %in% namespace$entities[[modelname]]
    if (isLocalEntity) {
		return(omxIdentifier(modelname, identifier))
    } else {
		return(identifier)
	}
}

omxConvertLabel <- function(label, modelname, dataname, namespace) {
	components <- unlist(strsplit(label, omxSeparatorChar, fixed = TRUE))
	if (omxIsDefinitionVariable(label)) {
		if (length(components) == 3) {
			return(label)
		} else {
			datasource <- unlist(strsplit(dataname, omxSeparatorChar, fixed = TRUE))[[1]]
			return(omxIdentifier(datasource, label))
		}
	}
	return(omxConvertIdentifier(label, modelname, namespace))
}

namespaceConvertMatrix <- function(matrix, modelname, dataname, namespace) {
	matrix@name <- omxIdentifier(modelname, matrix@name)
	free <- matrix@free
	select <- free == FALSE
	if (any(select)) {
		rows <- row(free)[select]
		cols <- col(free)[select]
		refNames <- matrix@labels[select]
		for(j in 1:length(refNames)) {
			row <- rows[j]
			col <- cols[j]
			matrix@labels[row,col] <- omxConvertLabel(refNames[j], modelname, dataname, namespace)
		}
	}
	return(matrix)
}

namespaceConvertAlgebra <- function(algebra, modelname, namespace) {
	algebra@name <- omxIdentifier(modelname, algebra@name)
	algebra@formula <- namespaceConvertFormula(algebra@formula, modelname, namespace)
	return(algebra)
}

namespaceConvertFormula <- function(formula, modelname, namespace) {
	if (length(formula) == 1) {
        if (is.symbol(formula) && 
            (as.character(formula) %in% namespace$parameters || 
             as.character(formula) %in% namespace$values)) {
        } else if (is.numeric(formula)) {
        } else {
            result <- omxConvertIdentifier(formula, modelname, namespace)
            if (is.symbol(formula) && is.character(result)) {
                formula <- as.symbol(result)
            } else {
                formula <- result
            }
        }
	} else {
		for (i in 2:length(formula)) {
			formula[[i]] <- namespaceConvertFormula(formula[[i]], modelname, namespace)
		}
	}
	return(formula)
}

namespaceConvertConstraint <- function(constraint, modelname, namespace) {
	constraint@name <- omxIdentifier(modelname, constraint@name)
	constraint@alg1 <- omxConvertIdentifier(constraint@alg1, modelname, namespace)
	constraint@alg2 <- omxConvertIdentifier(constraint@alg2, modelname, namespace)
	return(constraint)
}

namespaceConvertObjective <- function(objective, modelname, namespace) {
	if (!is.null(objective)) {
		objective <- omxObjFunNamespace(objective, modelname, namespace)
	}
	return(objective)
}

namespaceConvertData <- function(data, modelname) {
	if (!is.null(data)) {
		data@name <- omxIdentifier(modelname, data@name)
	}
	return(data)
}

namespaceSearch <- function(model, namespace, name, flat = FALSE) {
	if (namespace == model@name) {
		if (name %in% names(omxReservedNames) && 
			!is.null(omxReservedNames[[name]]@search)) {
			return(omxReservedNames[[name]]@search(model))
		}
		first <- model@matrices[[name]]
		second <- model@algebras[[name]]
		third <- model@submodels[[name]]
		fourth <- model@constraints[[name]]
		if (flat) { fifth <- model@objectives[[name]] }
		else { fifth <- NULL }
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
		} else {
            return(fifth)
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
		if (name %in% names(omxReservedNames) && 
			!is.null(omxReservedNames[[name]]@replace)) {
			return(omxReservedNames[[name]]@replace(model, value))
		}
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

