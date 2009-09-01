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
		options = "list",
		constMatrices = "list",
		output = "list"
))

omxModelTypes[['raw']] <- "MxModel"

setMethod("initialize", "MxModel",
	function(.Object, name = character(),
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
		.Object@options <- list()
		.Object@output <- list()
		.Object@constMatrices <- list()
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

setGeneric("omxVerifyModel", function(model) {
    return(standardGeneric("omxVerifyModel"))
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

setMethod("$", "MxModel",
	function(x, name) {
		return(omxExtractMethod(x, name))
	}
)

setReplaceMethod("$", "MxModel",
	function(x, name, value) {
		return(omxReplaceMethod(x, name, value))
	}
)


omxExtractMethod <- function(model, index) {
	pair <- omxReverseIdentifier(model, index)
	namespace <- pair[[1]]
	name <- pair[[2]]
	return(namespaceSearch(model, namespace, name))
}

omxReplaceMethod <- function(model, index, value) {
	pair <- omxReverseIdentifier(model, index)
	namespace <- pair[[1]]
	name <- pair[[2]]
	return(namespaceSearchReplace(model, namespace, name, value))
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
	retval <- firstArgument(model, name)
	first <- retval[[1]]
	model <- retval[[2]]
	name  <- retval[[3]]
	model <- typeArgument(model, type)
	lst <- c(first, list(...))
	lst <- mappendHelper(lst, list())
	model <- omxModelBuilder(model, lst, name, manifestVars,
		latentVars, remove, independent)
	model@options <- getOption('mxOptimizerOptions')
	return(model)
}

firstArgument <- function(model, name) {
	first <- NULL
	defaultType <- omxModelTypes[[getOption("mxDefaultType")]]
	if (is(model, "MxModel")) {
	} else {
		if (single.na(model)) {
		} else if (typeof(model) == "character") {
			name <- model
		} else if (isS4(model)) {
			first <- model
		} else {
			first <- list(model)
		}
		if (is.na(name)) {
			name <- omxUntitledName()
		}
		omxVerifyName(name)
		model <- new(defaultType, name = name)
	}
	return(list(first, model, name))
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
		model <- modelRemoveVariables(model, latentVars, manifestVars)
	} else if (length(manifestVars) + length(latentVars) > 0) {
		latentVars <- as.character(latentVars)
		manifestVars <- as.character(manifestVars)
		checkVariables(model, latentVars, manifestVars)
		model <- modelAddVariables(model, latentVars, manifestVars)
	}
	return(model)
}

listArgument <- function(model, lst, remove) {
	if(remove == TRUE) {
		model <- modelRemoveEntries(model, lst)
	} else {
		model <- modelAddEntries(model, lst)
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

checkVariables <- function(model, latentVars, manifestVars) {
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
	if (length(unique(latentVars)) != length(latentVars)) {
		stop("The latent variables list contains duplicate elements",
			call. = FALSE)
	}
	if (length(unique(manifestVars)) != length(manifestVars)) {
		stop("The manifest variables list contains duplicate elements",
			call. = FALSE)
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

setMethod("omxVerifyModel", "MxModel", function(model) {
    return(TRUE)
})

# End implementation of generics

modelAddVariables <- function(model, latent, manifest) {
	model@latentVars   <- union(model@latentVars, latent)
	model@manifestVars <- union(model@manifestVars, manifest)
	return(model)
}

modelRemoveVariables <- function(model, latent, manifest) {
	model@latentVars <- setdiff(model@latentVars, latent)
	model@manifestVars <- setdiff(model@manifestVars, manifest)
	return(model)
}
	
modelAddEntries <- function(model, entries) {
	if (length(entries) == 0) {
		return(model)
	}
	tuple <- modelAddFilter(model, entries, list(), list())
	namedEntities <- tuple[[1]]
	bounds        <- tuple[[2]]
	if (length(namedEntities) > 0) for(i in 1:length(namedEntities)) {
		model <- addSingleNamedEntity(model, namedEntities[[i]])
	}
	model <- modelAddBounds(model, bounds)
	return(model)
}

modelRemoveEntries <- function(model, entries) {
	if (length(entries) == 0) {
		return(model)
	}
	tuple <- modelRemoveFilter(model, entries, list(), list())
	namedEntities <- tuple[[1]]
	bounds        <- tuple[[2]]
	if (length(namedEntities) > 0) for(i in 1:length(namedEntities)) {
		model <- removeSingleNamedEntity(model, namedEntities[[i]])
	}
	model <- modelRemoveBounds(model, bounds)
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

modelAddFilter <- function(model, entries, namedEntities, bounds) {
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
	return(modelAddFilter(model, entries[-1], namedEntities, bounds))
}

modelRemoveFilter <- function(model, entries, names, bounds) {
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
	return(modelRemoveFilter(model, entries[-1], names, bounds))
}

addSingleNamedEntity <- function(model, entity) {
	model[[entity@name]] <- entity
	return(model)
}

removeSingleNamedEntity <- function(model, name) {
	model[[name]] <- NULL
	return(model)
}
