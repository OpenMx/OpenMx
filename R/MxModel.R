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


setClass(Class = "MxModel",
	representation = representation(
		name = "character",
		matrices = "list",
		algebras = "list",
		constraints = "list",
		latentVars = "character",
		manifestVars = "character",
		data = "MxData",
		submodels = "list",
		objective = "MxObjective",
		independent = "logical",
		output = "list"
))

omxModelTypes[['raw']] <- "MxModel"

setMethod("initialize", "MxModel",
	function(.Object, name = omxUntitledName(),  
		latentVars = character(), manifestVars = character(), 
		matrices = list(), algebras = list(), 
		constraints = list(), data = NULL, submodels = list(), 
		objective = NULL, independent = FALSE) {
		.Object@name <- name
		.Object@latentVars <- latentVars
		.Object@manifestVars <- manifestVars
		.Object@matrices <- matrices
		.Object@algebras <- algebras
		.Object@constraints <- constraints
		.Object@data <- data
		.Object@submodels <- submodels
		.Object@objective <- objective
		.Object@independent <- independent
		.Object@output <- list()
		.Object <- omxInitModel(.Object)
		return(.Object)
	}
)

# Begin declaration of generics

setGeneric("omxInitModel", function(model) {
	return(standardGeneric("omxInitModel")) } )

setGeneric("omxModelBuilder", function(model, lst, name, 
	manifestVars, latentVars, remove, independent) {
	return(standardGeneric("omxModelBuilder")) } )

setGeneric("omxTypeName", function(model) { 
	return(standardGeneric("omxTypeName")) 
})

# End declaration of generics

setMethod("[[", "MxModel",
	function(x, i, j, ..., drop = FALSE) {
		return(omxExtractMethod(x, i))
	}
)

setReplaceMethod("[[", "MxModel",
	function(x, i, j, value) {
		return(omxReplaceMethod(x, i, value))
	}
)

omxExtractMethod <- function(model, index) {
	first <- model@matrices[[index]]
	second <- model@algebras[[index]]
	third <- model@submodels[[index]]
	fourth <- model@constraints[[index]]
	if (!is.null(model@objective) && index == model@objective@name) {
		return(model@objective)
	} else if (!is.null(model@data) && index == model@data@name) {
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
}

omxReplaceMethod <- function(model, index, value) {
	current <- model[[index]]
	if (is.null(current) && is.null(value)) {
		return(model)
	}
	if(index == model@name) {
		stop(paste(omxQuotes(index), 
			"is already used as the name of the model"), call. = FALSE)
	}
	if(!is.null(current) && !is.null(value) && 
			!omxSameType(current, value)) {
		stop(paste("There already exists an object", 
				omxQuotes(index), 
				"in this model of different type"), call. = FALSE)
	}
	if(!is.null(value)) {
		value@name <- index
		test <- value
	} else {
		test <- current
	}
	if (is(test,"MxMatrix")) {
		model@matrices[[index]] <- value
	} else if (is(test,"MxAlgebra")) {
		model@algebras[[index]] <- value
	} else if (is(test,"MxModel")) {
		model@submodels[[index]] <- value	
	} else if (is(test,"MxObjective")) {
		model@objective <- value
	} else if (is(test,"MxData")) {
		model@data <- value
	} else if (is(test,"MxConstraint")) {
		model@constraints[[index]] <- value
	} else {
		stop(paste(test, "is of unknown value for replacement using name",
			index, "in model", model@name), call. = FALSE)
	}
	return(model)
}

omxSameType <- function(a, b) {
	return( (is(a, "MxModel") && is(b, "MxModel")) ||
			(is(a, "MxMatrix") && is(b, "MxMatrix")) ||
			(is(a, "MxAlgebra") && is(b, "MxAlgebra")) ||
			(is(a, "MxObjective") && is(b, "MxObjective")) ||
			(is(a, "MxConstraint") && is(b, "MxConstraint")) ||
			(is(a, "MxData") && is(b, "MxData")))
}

mxModel <- function(model = NA, ..., manifestVars = NA, latentVars = NA,
	remove = FALSE, independent = NA, type = NA, name = NA) {
	retval <- firstArgument(model)
	first <- retval[[1]]
	model <- retval[[2]]
	model <- typeArgument(model, type)
	lst <- c(first, list(...))
	lst <- mappendHelper(lst, list())
	model <- omxModelBuilder(model, lst, name, manifestVars,
		latentVars, remove, independent)
	return(model)
}
firstArgument <- function(model) {
	first <- NULL
	defaultType <- omxModelTypes[[getOption("mxDefaultType")]]
	if (typeof(model) != "S4" && is.na(model)) {
		model <- new(defaultType)	
	} else if (typeof(model) == "character") {
		model <- new(defaultType, name = model)
	} else if(!is(model, "MxModel")) {
		if(isS4(model)) {
	 		first <- model
		} else {
			first <- list(model)
		}
		model <- new(defaultType)
	}
	return(list(first, model))
}

typeArgument <- function(model, type) {
	if (!is.na(type)) {
		if (is.null(omxModelTypes[[type]])) {
			stop(paste("The model type", omxQuotes(type), 
				"is not in the the list of acceptable types:",
				omxQuotes(names(omxModelTypes))), call. = FALSE)
		}
		typename <- omxModelTypes[[type]]
		class(model) <- typename
		model <- omxInitModel(model)
	}
	return(model)
}

omxGenericModelBuilder <- function(model, lst, name, 
	manifestVars, latentVars, remove, independent) {
	model <- variablesArgument(model, manifestVars, latentVars, remove)
	model <- listArgument(model, lst, remove)
	model <- independentArgument(model, independent)
	model <- nameArgument(model, name)
	return(model)
}

variablesArgument <- function(model, manifestVars, latentVars, remove) {
	if (single.na(manifestVars)) {
		manifestVars <- character()
	}
	if (single.na(latentVars)) {
		latentVars <- character()
	}
	if (remove == TRUE) {

	} else if (length(manifestVars) + length(latentVars) > 0) {
		latentVars <- as.character(latentVars)
		manifestVars <- as.character(manifestVars)
		omxCheckVariables(model, latentVars, manifestVars)
		model <- omxAddVariables(model, latentVars, manifestVars)
	}
	return(model)
}

listArgument <- function(model, lst, remove) {
	if(remove == TRUE) {
		model <- omxRemoveEntries(model, lst)
	} else {
		model <- omxAddEntries(model, lst)
	}
	return(model)
}

independentArgument <- function(model, independent) {
	if(!is.na(independent)) {
		model@independent <- independent
	}
	return(model)
}

nameArgument <- function(model, name) {
	if(!is.na(name)) {
		model@name <- name
	}
	return(model)
}

omxCheckVariables <- function(model, latentVars, manifestVars) {
	common <- intersect(latentVars, manifestVars)
	if (length(common) > 0) {
		stop(paste("The following variables cannot",
			"be both latent and manifest:",
			omxQuotes(common)), call. = FALSE)
	}
	common <- intersect(model@latentVars, manifestVars)
	if (length(common) > 0) {
		stop(paste("The following variables cannot",
			"be both latent and manifest:",
			omxQuotes(common)), call. = FALSE)
	}
	common <- intersect(model@manifestVars, latentVars)
	if (length(common) > 0) {
		stop(paste("The following variables cannot",
			"be both latent and manifest",
			omxQuotes(common)), call. = FALSE)
	}
	if (any(is.na(latentVars))) {
		stop("NA is not allowed as a latent variable", call. = FALSE)
	}
	if (any(is.na(manifestVars))) {
		stop("NA is not allowed as a manifest variable", call. = FALSE)
	}
}

# Begin implementation of generics

setMethod("omxModelBuilder", "MxModel", omxGenericModelBuilder)

setMethod("omxInitModel", "MxModel", function(model) { 
	return(model)
})

setMethod("omxTypeName", "MxModel", function(model) { 
	return("unspecified")
})

# End implementation of generics

omxAddVariables <- function(model, latent, manifest) {
	model@latentVars   <- union(model@latentVars, latent)
	model@manifestVars <- union(model@manifestVars, manifest)
	return(model)
}
	
omxAddEntries <- function(model, entries) {
	if (length(entries) == 0) {
		return(model)
	}
	tuple <- omxModelAddFilter(model, entries, list(), list())
	namedEntities <- tuple[[1]]
	bounds        <- tuple[[2]]
	if (length(namedEntities) > 0) for(i in 1:length(namedEntities)) {
		model <- omxAddSingleNamedEntity(model, namedEntities[[i]])
	}
	model <- omxAddBounds(model, bounds)
	return(model)
}

omxRemoveEntries <- function(model, entries) {
	if (length(entries) == 0) {
		return(model)
	}
	tuple <- omxModelRemoveFilter(model, entries, list(), list())
	namedEntities <- tuple[[1]]
	bounds        <- tuple[[2]]
	if (length(namedEntities) > 0) for(i in 1:length(namedEntities)) {
		model <- omxRemoveSingleNamedEntity(model, namedEntities[[i]])
	}
	model <- omxRemoveBounds(model, bounds)
	return(model)
}

mappendHelper <- function(lst, result) {
	if (length(lst) == 0) {
		return(result)
	} else if (length(lst) == 1) {
		len <- length(result)
		addition <- lst[[1]]
		if (is.list(addition)) {
			result <- append(result, addition)
		} else {
			result[[len + 1]] <- addition
		}
		return(result)
	} else {
		len <- length(result)
		addition <- lst[[1]]
		if (is.list(addition)) {
			result <- append(result, addition)
		} else {
			result[[len + 1]] <- addition
		}
		return(mappendHelper(lst[2:length(lst)], result))
	}
}

omxModelAddFilter <- function(model, entries, namedEntities, bounds) {
	if (length(entries) == 0) {
		return(list(namedEntities, bounds))
	}
	head <- entries[[1]]
	nLength <- length(namedEntities)
	bLength <- length(bounds)
	if (is.null(head)) {
	} else if(isS4(head) && ("name" %in% slotNames(head))) {
		namedEntities[[nLength + 1]] <- head
	} else if(is(head, "MxBounds")) {
		bounds[[bLength + 1]] <- head
	} else if(omxIsPath(head)) {
		stop(paste("The model type of model",
			omxQuotes(model@name), "does not recognize paths."),
			call. = FALSE)
	} else {
		stop(paste("Cannot add the following item into the model:", 
			head), call. = FALSE)
	}
	return(omxModelAddFilter(model, entries[-1], namedEntities, bounds))
}

omxModelRemoveFilter <- function(model, entries, names, bounds) {
	if (length(entries) == 0) {
		return(list(names, bounds))
	}
	head <- entries[[1]]
	nLength <- length(names)
	bLength <- length(bounds)
	if (is.null(head)) {
	} else if(is.character(head) && (length(head) == 1)) {
		names[[nLength + 1]] <- head
	} else if(is(head, "MxBounds")) {
		bounds[[bLength + 1]] <- head
	} else if(omxIsPath(head)) {
		stop(paste("The model type of model",
			omxQuotes(model@name), "does not recognize paths."),
			call. = FALSE)
	} else {
		stop(paste("Cannot remove the following item from the model:", 
			head), call. = FALSE)
	}
	return(omxModelRemoveFilter(model, entries[-1], names, bounds))
}

omxAddSingleNamedEntity <- function(model, entity) {
	model[[entity@name]] <- entity
	return(model)
}

omxRemoveSingleNamedEntity <- function(model, name) {
	model[[name]] <- NULL
	return(model)
}
