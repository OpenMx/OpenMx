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


setClass(Class = "MxMLObjective",
	representation = representation(
		covariance = "MxCharOrNumber",
		means = "MxCharOrNumber",
		dims = "character",
		thresholds = "character"),
	contains = "MxBaseObjective")

setMethod("initialize", "MxMLObjective",
	function(.Object, covariance, means, dims, thresholds, 
		data = as.integer(NA), name = 'objective') {
		.Object@name <- name
		.Object@covariance <- covariance
		.Object@means <- means
		.Object@data <- data
		.Object@dims <- dims
		.Object@thresholds <- thresholds
		return(.Object)
	}
)

setMethod("genericObjDependencies", signature("MxMLObjective"),
	function(.Object, dependencies) {
	sources <- c(.Object@covariance, .Object@means, .Object@thresholds)
	sources <- sources[!is.na(sources)]
	dependencies <- omxAddDependency(sources, .Object@name, dependencies)
	return(dependencies)
})


setMethod("genericObjFunNamespace", signature("MxMLObjective"), 
	function(.Object, modelname, namespace) {
		.Object@name <- omxIdentifier(modelname, .Object@name)
		.Object@covariance <- omxConvertIdentifier(.Object@covariance, 
			modelname, namespace)
		.Object@means <- omxConvertIdentifier(.Object@means, 
			modelname, namespace)
		.Object@data <- omxConvertIdentifier(.Object@data, 
			modelname, namespace)
		.Object@thresholds <- sapply(.Object@thresholds, 
			omxConvertIdentifier, modelname, namespace)			
		return(.Object)
})

setMethod("genericObjRename", signature("MxMLObjective"),
	function(.Object, oldname, newname) {
		.Object@covariance <- renameReference(.Object@covariance, oldname, newname)
		.Object@means <- renameReference(.Object@means, oldname, newname)		
		.Object@data <- renameReference(.Object@data, oldname, newname)	
		.Object@thresholds <- sapply(.Object@thresholds, renameReference, oldname, newname)		
		return(.Object)
})

setMethod("genericObjFunConvert", signature("MxMLObjective"), 
	function(.Object, flatModel, model, defVars) {
		modelname <- omxReverseIdentifier(model, .Object@name)[[1]]
		name <- .Object@name
		covariance <- .Object@covariance
		means <- .Object@means
		data <- .Object@data
		dims <- .Object@dims
		thresholds <- .Object@thresholds
		covarianceIndex <- omxLocateIndex(flatModel, covariance, name)
		if(is.na(data)) {
			msg <- paste("In model", omxQuotes(modelname),
				"the ML objective does not have a dataset specified")
			stop(msg, call. = FALSE)
		}
		mxDataObject <- flatModel@datasets[[data]]
		verifyObservedNames(mxDataObject@observed, mxDataObject@type, flatModel, modelname, "ML")
		checkNumericData(mxDataObject)
		verifyExpectedNames(covariance, means, flatModel, modelname, "ML")
		meansIndex <- omxLocateIndex(flatModel, means, name)
		dIndex <- omxLocateIndex(flatModel, data, name)
		.Object@covariance <- covarianceIndex
		.Object@means <- meansIndex
		.Object@data <- dIndex
		if (single.na(.Object@dims)) {
			.Object@dims <- rownames(flatModel[[covariance]])
		}
		return(.Object)
})


setMethod("genericObjModelConvert", "MxMLObjective",
	function(.Object, job, model, flatJob) {
		if(is.na(.Object@data)) {
			msg <- paste("The ML objective",
				"does not have a dataset associated with it in model",
				omxQuotes(model@name))
			stop(msg, call. = FALSE)
		}
		job <- updateObjectiveDimnames(.Object, job, model@name, "ML")
		if (flatJob@datasets[[.Object@data]]@type != 'raw') {
			return(job)
		}
		if (is.na(.Object@means)) {
			msg <- paste("In model", omxQuotes(model@name),
				"the ML objective has a raw dataset specified",
				"but no expected means vector")
			stop(msg, call. = FALSE)
		}
		objective <- mxFIMLObjective(.Object@covariance, .Object@means, .Object@thresholds)
		job[[model@name]][['objective']] <- objective
		return(job)
	}
)


mxMLObjective <- function(covariance, means = NA, dimnames = NA, thresholds = NA) {
	if (missing(covariance) || typeof(covariance) != "character") {
		stop("Covariance argument is not a string (the name of the expected covariance matrix)")
	}
	if (!(single.na(means) || typeof(means) == "character")) {
		stop("Means argument is not a string (the name of the expected means matrix)")
	}
	if (is.na(means)) means <- as.integer(NA)
	if (single.na(thresholds)) thresholds <- as.character(NA)
	if (single.na(dimnames)) dimnames <- as.character(NA)
	if (!is.vector(dimnames) || typeof(dimnames) != 'character') {
		stop("Dimnames argument is not a character vector")
	}
	if (length(thresholds) != 1) {
		stop("Thresholds argument must be a single matrix or algebra name")
	}
	if (length(dimnames) == 0) {
		stop("Dimnames argument cannot be an empty vector")
	}
	if (length(dimnames) > 1 && any(is.na(dimnames))) {
		stop("NA values are not allowed for dimnames vector")
	}
	return(new("MxMLObjective", covariance, means, dimnames, thresholds))
}

displayMLObjective <- function(objective) {
	cat("MxMLObjective", omxQuotes(objective@name), '\n')
	cat("@covariance :", omxQuotes(objective@covariance), '\n')
	if (single.na(objective@means)) {
		cat("@means : NA \n")
	} else {
		cat("@means :", omxQuotes(objective@means), '\n')
	}
	if (single.na(objective@thresholds)) {
		cat("@thresholds : NA \n")
	} else {
		cat("@thresholds :", omxQuotes(objective@thresholds), '\n')
	}
	if (single.na(objective@dims)) {
		cat("@dims : NA \n")
	} else {
		cat("@dims :", omxQuotes(objective@dims), '\n')
	}
	if (length(objective@result) == 0) {
		cat("@result: (not yet computed) ")
	} else {
		cat("@result:\n")
	}
	print(objective@result)
	invisible(objective)
}


setMethod("print", "MxMLObjective", function(x,...) { 
	displayMLObjective(x) 
})

setMethod("show", "MxMLObjective", function(object) { 
	displayMLObjective(object) 
})
