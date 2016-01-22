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

setOldClass('package_version')
setClass(Class = "MxModel",
	representation = representation(
		name = "character",
		matrices = "list",
		algebras = "list",
		constraints = "list",
		intervals = "list",
		latentVars = "MxCharOrList",
		manifestVars = "MxCharOrList",
		data = "MxData",
		submodels = "list",
		expectation = "MxExpectation",
		fitfunction = "MxFitFunction",
		compute = "MxCompute",
		independent = "logical",
		options = "list",
		output = "list",
		runstate = "list",
		.newobjects = "logical",
		.resetdata = "logical",
	        .wasRun = "logical",
	    .modifiedSinceRun = "logical",
	    .version = "package_version"
))

imxModelTypes[['default']] <- "MxModel"

setMethod("initialize", "MxModel",
	function(.Object, name = character()) {
		.Object@name <- name
		.Object@latentVars <- character()
		.Object@manifestVars <- character()
		.Object@matrices <- list()
		.Object@algebras <- list()
		.Object@constraints <- list()
		.Object@data <- NULL
		.Object@submodels <- list()
		.Object@expectation <- NULL
		.Object@fitfunction <- NULL
		.Object@compute <- NULL
		.Object@independent <- FALSE
		.Object@options <- list()
		.Object@output <- list()
		.Object@runstate <- list()
		.Object@.newobjects <- FALSE
		.Object@.resetdata <- FALSE
	        .Object@.wasRun <- FALSE
	        .Object@.modifiedSinceRun <- FALSE
		.Object@.version <- packageVersion("OpenMx")
		.Object <- imxInitModel(.Object)
		return(.Object)
	}
)

##' imxInitModel
##'
##' This is an internal function exported for those people who know
##' what they are doing.
##'
##' @param model model
##' @aliases
##' imxInitModel,MxModel-method
##' imxInitModel,MxRAMModel-method
##' imxInitModel,MxLISRELModel-method
setGeneric("imxInitModel", function(model) {
	return(standardGeneric("imxInitModel")) } )

##' imxModelBuilder
##'
##' This is an internal function exported for those people who know
##' what they are doing.
##' 
##' TODO: It probably makes sense to split this into separate
##' methods. For example, modelAddVariables and modelRemoveVariables
##' could be their own methods. This would reduce some cut&paste
##' duplication.
##'
##' @param model model
##' @param lst lst
##' @param name name
##' @param manifestVars manifestVars
##' @param latentVars latentVars
##' @param submodels submodels
##' @param remove remove
##' @param independent independent
##' @aliases
##' imxModelBuilder,MxLISRELModel-method
##' imxModelBuilder,MxModel-method
##' imxModelBuilder,MxRAMModel-method
setGeneric("imxModelBuilder", function(model, lst, name, 
	manifestVars, latentVars, submodels, remove, independent) {
	return(standardGeneric("imxModelBuilder")) } )

##' imxTypeName
##'
##' This is an internal function exported for those people who know
##' what they are doing.
##'
##' @param model model
##' @aliases
##' imxTypeName,MxLISRELModel-method
##' imxTypeName,MxModel-method
##' imxTypeName,MxRAMModel-method
setGeneric("imxTypeName", function(model) { 
	return(standardGeneric("imxTypeName")) 
})

##' imxVerifyModel
##' 
##' This is an internal function exported for those people who know
##' what they are doing.
##'
##' @param model model
##' @aliases
##' imxVerifyModel,MxLISRELModel-method
##' imxVerifyModel,MxModel-method
##' imxVerifyModel,MxRAMModel-method
setGeneric("imxVerifyModel", function(model) {
    return(standardGeneric("imxVerifyModel"))
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
	if (!is.null(model@fitfunction)) {
		retval <- union(retval, model@fitfunction@name)
	}
	if (!is.null(model@expectation)) {
		retval <- union(retval, model@expectation@name)
	}
	if (!is.null(model@data)) {
		retval <- union(retval, model@data@name)
	}
	return(retval)
}

setMethod("[[", "MxModel",
	function(x, i, j, ..., drop = FALSE) {
		return(imxExtractMethod(x, i))
	}
)

setReplaceMethod("[[", "MxModel",
	function(x, i, j, value) {
		return(imxReplaceMethod(x, i, value))
	}
)

# These are slots that are intended to be directly viewable by the user.
# Included separately so that they are the same between the $ and names() operators.
publicMxModelSlots <- c("name", "matrices", "algebras", "constraints", "data", "submodels", "output", "compute", "options", "intervals")
visibleMxModelSlots <- c("name", "options", "compute", "output", "intervals")

setMethod("$", "MxModel",
	function(x, name) {
        result <- imxExtractMethod(x, name)
        if(name %in% publicMxModelSlots) {
            result <- imxExtractSlot(x, name)
        }
		return(result)
	}
)

setMethod("names", "MxModel",
	function(x) {
		submodels <- names(x@submodels)
		locals <- generateLocalNames(x)
		slots <- publicMxModelSlots
		output <- c(submodels, locals, slots)
		output <- gsub("(\\w+\\W+.*)", "'\\1'", output)
		return(output)
	}
)

.DollarNames.MxModel <- function(x, pattern) {   # .DollarNames is an S3 Generic
		submodels <- names(x@submodels)
		locals <- generateLocalNames(x)
		slots <- imxGetSlotDisplayNames(x, slotList=visibleMxModelSlots)
		output <- c(submodels, locals, slots)
		output <- gsub("(\\w+\\W+.*)", "'\\1'", output)
		return(grep(pattern, output, value=TRUE))
	}

setReplaceMethod("$", "MxModel",
	function(x, name, value) {
		if(name == "output") {
			stop("You cannot directly set the output of a model.  Use mxRun() if you want output.")
		} 
		if(name == "name") {
			stop("You cannot directly set the name of a model.  To rename the model, use model<-mxModel(model, name=\"NewName\").")
		}
		if(name %in% c("matrices", "algebras", "submodels")) {
			stop(paste("You cannot directly set the", name,
				   "of a model.  To set objects in the model, use the mxModel() function."))
		}
		return(imxReplaceMethod(x, name, value))
	}
)


##' imxSameType
##' 
##' This is an internal function exported for those people who know
##' what they are doing.
##'
##' @param a a
##' @param b b
imxSameType <- function(a, b) {
	return( (is(a, "MxModel") && is(b, "MxModel")) ||
			(is(a, "MxMatrix") && is(b, "MxMatrix")) ||
			(is(a, "MxAlgebra") && is(b, "MxAlgebra")) ||
			(is(a, "MxExpectation") && is(b, "MxExpectation")) ||
			(is(a, "MxFitFunction") && is(b, "MxFitFunction")) ||
			(is(a, "MxCompute") && is(b, "MxCompute")) ||
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
	filter <- sapply(lst, is, "MxModel")
	submodels <- lst[filter]
	lst <- lst[!filter]
	model <- imxModelBuilder(model, lst, name, manifestVars,
		latentVars, submodels, remove, independent)
	return(model)
}

firstArgument <- function(model, name) {
	first <- NULL
	defaultType <- imxModelTypes[[getOption("mxDefaultType")]]
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
			name <- imxUntitledName()
		}
		imxVerifyName(name, -1)
		model <- new(defaultType, name)
	}
	return(list(first, model, name))
}

typeArgument <- function(model, type) {
	if (!is.na(type)) {
		if (is.null(imxModelTypes[[type]])) {
			stop(paste("The model type", omxQuotes(type), 
				"is not in the the list of acceptable types:",
				omxQuotes(names(imxModelTypes))), call. = FALSE)
		}
		typename <- imxModelTypes[[type]]
		class(model) <- typename
		model <- imxInitModel(model)
	}
	return(model)
}

##' imxGenericModelBuilder
##'
##' This is an internal function exported for those people who know
##' what they are doing.
##'
##' @param model model
##' @param lst lst
##' @param name name
##' @param manifestVars manifestVars
##' @param latentVars latentVars
##' @param submodels submodels
##' @param remove remove
##' @param independent independent
imxGenericModelBuilder <- function(model, lst, name, 
	manifestVars, latentVars, submodels, remove, independent) {
	model <- nameArgument(model, name)
	model <- variablesArgument(model, manifestVars, latentVars, submodels, remove)
	model <- listArgument(model, lst, remove)
	model <- independentArgument(model, independent)
	return(model)
}

varsToCharacter <- function(vars, vartype) {
	if (is.list(vars)) {
		varnames <- names(vars)
		if (length(varnames) == 0) {
			return(as.character(unlist(vars)))	
		} else {
			result <- pmatch(varnames, imxVariableTypes)
			illegal <- which(is.na(result))
			if (length(illegal) > 0) {
				if (length(illegal) == 1) {
					ctgMsg <- "category"
				} else {
					ctgMsg <- "categories"
				}
				msg <- paste("In the", vartype, "variables",
					"the", ctgMsg,
					omxQuotes(varnames[illegal]), "did not match",
					"to a valid category or two categories matched",
					"to the same string (see 'imxVariableTypes'",
					"for the list of legal categories)")
				stop(msg, call. = FALSE)
			}
			varnames <- imxVariableTypes[result]
			vars <- lapply(vars, as.character)
			names(vars) <- varnames
			return(vars)
		}
	} else {
		return(as.character(vars))
	}
}

variablesArgument <- function(model, manifestVars, latentVars, submodels, remove) {
	if (single.na(manifestVars)) {
		manifestVars <- character()
	}
	if (single.na(latentVars)) {
		latentVars <- character()
	}
	if (remove == TRUE) {
		model <- modelRemoveVariables(model, latentVars, manifestVars)
		if (length(submodels)) for(i in 1:length(submodels)) {
			model <- removeSingleNamedEntity(model, submodels[[i]])
		}
	} else {
		if (length(manifestVars) + length(latentVars) > 0) {
			latentVars <- varsToCharacter(latentVars, "latent")
			manifestVars <- varsToCharacter(manifestVars, "manifest")
			checkVariables(model, latentVars, manifestVars)
			model <- modelAddVariables(model, latentVars, manifestVars)
		}
		if (length(submodels)) for(i in 1:length(submodels)) {
			model <- addSingleNamedEntity(model, submodels[[i]])
		}
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

checkVariables <- function(model, latentVars, manifestVars, submodels) {
	latentVars   <- unlist(latentVars, use.names = FALSE)
	manifestVars <- unlist(manifestVars, use.names = FALSE)
	modelLatent <- unlist(model@latentVars, use.names = FALSE)
	modelManifest <- unlist(model@manifestVars, use.names = FALSE)
	common <- intersect(latentVars, manifestVars)    
	if (length(common) > 0) {
		stop(paste("The following variables cannot",
			"be both latent and manifest:",
			omxQuotes(common)), call. = FALSE)
	}
	common <- intersect(modelLatent, manifestVars)
	if (length(common) > 0) {
		stop(paste("The following variables cannot",
			"be both latent and manifest:",
			omxQuotes(common)), call. = FALSE)
	}
	common <- intersect(modelManifest, latentVars)
	if (length(common) > 0) {
		stop(paste("The following variables cannot",
			"be both latent and manifest",
			omxQuotes(common)), call. = FALSE)
	}
	common <- intersect(modelManifest, manifestVars)
	if (length(common) > 0) {
		stop(paste("The following manifest variables",
			"have already been declared",
			omxQuotes(common)), call. = FALSE)
	}
	common <- intersect(modelLatent, latentVars)
	if (length(common) > 0) {
		stop(paste("The following latent variables",
			"have already been declared",
			omxQuotes(common)), call. = FALSE)
	}
	if (any(is.na(latentVars))) {
		stop("NA is not allowed as a latent variable", call. = FALSE)
	}
	if (any(is.na(manifestVars))) {
		stop("NA is not allowed as a manifest variable", call. = FALSE)
	}
	if (length(unique(latentVars)) != length(latentVars)) {
		stop(paste("The following variables in the latentVars list are duplicated:", 
		omxQuotes(latentVars[duplicated(latentVars)])), call. = FALSE)
	}
	if (length(unique(manifestVars)) != length(manifestVars)) {
		stop(paste("The following variables in the manifestVars list are duplicated:", 
		omxQuotes(manifestVars[duplicated(manifestVars)])), call. = FALSE)
	}
}

# Begin implementation of generics

setMethod("imxModelBuilder", "MxModel", imxGenericModelBuilder)

setMethod("imxInitModel", "MxModel", function(model) { 
	return(model)
})

setMethod("imxTypeName", "MxModel", function(model) { 
	return("default")
})

setMethod("imxVerifyModel", "MxModel", function(model) {
    return(TRUE)
})

# End implementation of generics

addVariablesHelper <- function(model, vartype, vars) {
	modelvars <- slot(model, vartype)

	if (length(vars) == 0) {
		return(model)
	} else if (length(modelvars) == 0) {
		slot(model, vartype) <- vars
		return(model)
	}

	if (is.list(vars) && !is.list(modelvars)) {
		msg <- paste("The", vartype, "variables in",
			"the call to mxModel() have been separated",
			"into categories, and the existing", vartype,
			"variables do not have categories.")
		stop(msg, call. = FALSE)
	} else if (!is.list(vars) && is.list(modelvars)) {
		msg <- paste("The", vartype, "variables in",
			"the call to mxModel() have not been separated",
			"into categories, and the existing", vartype,
			"variables do have categories.")
		stop(msg, call. = FALSE)
	}

	if (is.character(vars) && is.character(modelvars)) {
		modelvars <- c(modelvars, vars)
		slot(model, vartype) <- modelvars
	} else {
		varnames <- names(vars)
		offsets <- pmatch(varnames, imxVariableTypes)
		for(i in 1:length(vars)) {
			currOffset <- offsets[[i]]
			currName <- imxVariableTypes[[currOffset]]
			modelvars[[currName]] <- c(modelvars[[currName]], vars[[i]])
		}
		slot(model, vartype) <- modelvars
	}

	return(model)
}

modelAddVariables <- function(model, latent, manifest) {
	model <- addVariablesHelper(model, "latentVars", latent)
	model <- addVariablesHelper(model, "manifestVars", manifest)
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
	intervals     <- expandIntervals(intervals)
	names(intervals) <- sapply(intervals, slot, "reference")
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
	intervals     <- expandIntervals(intervals)	
	names(intervals) <- sapply(intervals, slot, "reference")
	if (length(namedEntities) > 0) for(i in 1:length(namedEntities)) {
		model <- removeSingleNamedEntity(model, namedEntities[[i]])
	}
	model <- modelRemoveBounds(model, bounds)
	model <- modelRemoveIntervals(model, intervals)
	return(model)
}

actionCorrespondingPredicate <- c('add' = 'into', 'remove' = 'from')

modelModifyFilter <- function(model, entries, action) {
	boundsFilter <- sapply(entries, is, "MxBounds")
	intervalFilter <- sapply(entries, is, "MxInterval")
	namedEntityFilter <- sapply(entries, function(x) {"name" %in% slotNames(x)})
	characterFilter <- sapply(entries, is.character)
	pathFilter <- sapply(entries, is, "MxPath")
  thresholdFilter <- sapply(entries, is, "MxThreshold")
	unknownFilter <- !(boundsFilter | namedEntityFilter | intervalFilter | characterFilter | thresholdFilter)
	if (any(pathFilter)) {
		stop(paste("The model type of model",
			omxQuotes(model@name), "does not recognize paths."),
			call. = FALSE)
	}
	if (any(thresholdFilter)) {
	  stop(paste("The model type of model",
	             omxQuotes(model@name), "does not recognize thresholds."),
	       call. = FALSE)
	}
	if (any(unknownFilter)) {
		stop(paste("Cannot", action, "the following item(s)", 
			actionCorrespondingPredicate[[action]], "the model:", 
			omxQuotes(sapply(entries[unknownFilter], deparse))), call. = FALSE)
	}
	if (any(namedEntityFilter) && action == 'remove') {
		stop(paste("Cannot use named entities when remove = TRUE.",
			"Instead give the name of the entity when removing it.",
			"See http://openmx.psyc.virginia.edu/wiki/mxmodel-help#Remove_an_object_from_a_model"))
	}
	if (any(characterFilter) && action == 'add') {
		stop(paste("I don't know what to do with the following strings",
			omxQuotes(entries[characterFilter]),
			"that have been passed into the function:",
			deparse(width.cutoff = 400L, imxLocateFunction("mxModel"))), call. = FALSE)
	}
	if (identical(action, 'add')) {
		return(list(entries[namedEntityFilter], entries[boundsFilter], entries[intervalFilter]))
	} else if (identical(action, 'remove')) {
		return(list(entries[characterFilter], entries[boundsFilter], entries[intervalFilter]))
	} else {
		stop(paste("Internal error, unidentified action:", omxQuotes(action)))
	}
}

addSingleNamedEntity <- function(model, entity) {
	if (!nzchar(entity@name)) {
		stop(paste("Entity",
			omxQuotes(class(entity)), "in model",
			omxQuotes(model@name), "needs a name"), call. = FALSE)
	}
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

setMethod("imxVerifyModel", "MxModel",
    function(model) {
        if (length(model@submodels) > 0) {
        	return(all(sapply(model@submodels, imxVerifyModel)))
        }
        return(TRUE)
    }
)
