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
		output = "list",
		runstate="list"
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
		.Object@runstate <- list()
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

generateParentNames <- function(model) {
	retval <- generateLocalNames(model)
	if (length(model@submodels) > 0) {
		retval <- union(retval, names(model@submodels))
		childNames <- unlist(lapply(model@submodels, generateChildNames))
		retval <- union(retval, childNames)
	}
	return(retval)
}

generateChildNames <- function(model) {
	retval <- generateLocalNames(model)	
	if (!is.null(retval)) {
		retval <- paste(model@name, retval, sep = ".")
	}
	if (length(model@submodels) > 0) {
		retval <- union(retval, names(model@submodels))
		childNames <- unlist(lapply(model@submodels, generateChildNames))
		retval <- union(retval, childNames)
	}
	return(retval)
}

generateLocalNames <- function(model) {
	matrices <- names(model@matrices)
	algebras <- names(model@algebras)
	constraints <- names(model@constraints)
	retval <- union(matrices, algebras)
	retval <- union(retval, constraints)
	if (!is.null(model@objective)) {
		retval <- union(retval, model@objective@name)
	}
	if (!is.null(model@data)) {
		retval <- union(retval, model@data@name)
	}
	return(retval)
}

setMethod("names", "MxModel",
	function(x) {
		generateParentNames(x)
	}
)

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
	return(namespaceSearch(model, index))
}

omxReplaceMethod <- function(model, index, value) {
	return(namespaceSearchReplace(model, index, value))
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
	lst <- unlist(lst)
	model <- omxModelBuilder(model, lst, name, manifestVars,
		latentVars, remove, independent)
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
		if (length(name) > 0 && is.na(name)) {
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
	return("default")
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
	tuple <- modelModifyFilter(model, entries, "add")
	namedEntities <- tuple[[1]]
	bounds        <- tuple[[2]]
	intervals     <- tuple[[3]]
	if (length(namedEntities) > 0) for(i in 1:length(namedEntities)) {
		model <- addSingleNamedEntity(model, namedEntities[[i]])
	}
	model <- modelAddBounds(model, bounds)
	model <- modelAddIntervals(model, intervals)
	return(model)
}

modelRemoveEntries <- function(model, entries) {
	if (length(entries) == 0) {
		return(model)
	}
	tuple <- modelModifyFilter(model, entries, "remove")
	namedEntities <- tuple[[1]]
	bounds        <- tuple[[2]]
	intervals     <- tuple[[3]]
	if (length(namedEntities) > 0) for(i in 1:length(namedEntities)) {
		model <- removeSingleNamedEntity(model, namedEntities[[i]])
	}
	model <- modelRemoveBounds(model, bounds)
	model <- modelRemoveIntervals(model, intervals)
	return(model)
}

modelModifyFilter <- function(model, entries, action) {
	boundsFilter <- sapply(entries, is, "MxBounds")
	intervalFilter <- sapply(entries, is, "MxInterval")
	namedFilter <- sapply(entries, function(x) {"name" %in% slotNames(x)})
	pathFilter <- sapply(entries, is, "MxPath")
	unknownFilter <- !(boundsFilter | namedFilter | intervalFilter)
	if (any(pathFilter)) {
		stop(paste("The model type of model",
			omxQuotes(model@name), "does not recognize paths."),
			call. = FALSE)
	}
	if (any(unknownFilter)) {
		stop(paste("Cannot", action, "the following item(s) into the model:", 
			entries[unknownFilter]), call. = FALSE)
	}
	return(list(entries[namedFilter], entries[boundsFilter], entries[intervalFilter]))
}

addSingleNamedEntity <- function(model, entity) {
	if (model@name == entity@name) {
		stop(paste("You cannot insert an entity named",
			omxQuotes(entity@name), "into a model named",
			omxQuotes(model@name)), call. = FALSE)
	}
	model[[entity@name]] <- entity
	return(model)
}

removeSingleNamedEntity <- function(model, name) {
	model[[name]] <- NULL
	return(model)
}
