#
#   Copyright 2007-2014 The OpenMx Project
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

setClass(Class = "MxExpectationNormal",
	representation = representation(
		covariance = "MxCharOrNumber",
		means = "MxCharOrNumber",
		definitionVars = "list",
		thresholds = "MxCharOrNumber",
		dims = "character",
		dataColumns = "numeric",
		thresholdColumns = "numeric",
		thresholdLevels = "numeric",
		threshnames = "character",
		ExpCov = "matrix",
		ExpMean = "matrix",
	        numStats = "numeric"),
	contains = "MxBaseExpectation")

setMethod("initialize", "MxExpectationNormal",
	function(.Object, covariance, means, dims, thresholds, threshnames,
		data = as.integer(NA), definitionVars = list(), name = 'expectation') {
		.Object@name <- name
		.Object@covariance <- covariance
		.Object@means <- means
		.Object@data <- data
		.Object@definitionVars <- definitionVars
		.Object@thresholds <- thresholds
		.Object@dims <- dims
		.Object@threshnames <- threshnames
		.Object@ExpCov <- matrix()
		.Object@ExpMean <- matrix()
		return(.Object)
	}
)

setMethod("qualifyNames", signature("MxExpectationNormal"), 
	function(.Object, modelname, namespace) {
		.Object@name <- imxIdentifier(modelname, .Object@name)
		.Object@covariance <- imxConvertIdentifier(.Object@covariance, 
			modelname, namespace)
		.Object@means <- imxConvertIdentifier(.Object@means, 
			modelname, namespace)
		.Object@data <- imxConvertIdentifier(.Object@data, 
			modelname, namespace)
		.Object@thresholds <- sapply(.Object@thresholds,
			imxConvertIdentifier, modelname, namespace)
		return(.Object)
})

setMethod("genericExpDependencies", signature("MxExpectationNormal"),
	function(.Object, dependencies) {
	sources <- c(.Object@covariance, .Object@means, .Object@thresholds)
	sources <- sources[!is.na(sources)]
	dependencies <- imxAddDependency(sources, .Object@name, dependencies)
	return(dependencies)
})

setMethod("genericExpAddEntities", "MxExpectationNormal",
	function(.Object, job, flatJob, labelsData) {
		precision <- "Function precision"
		if(!single.na(.Object@thresholds)) {
			if (is.null(job@options[[precision]])) {
				job <- mxOption(job, precision, 1e-9)
			}
		}
		return(job)
	}
)

setMethod("genericExpConvertEntities", "MxExpectationNormal",
	function(.Object, flatModel, namespace, labelsData) {
		flatModel <- updateExpectationDimnames(.Object, flatModel, labelsData)
		flatModel <- updateThresholdDimnames(.Object, flatModel, labelsData)
		return(flatModel)
	}
)

setMethod("genericExpRename", signature("MxExpectationNormal"),
	function(.Object, oldname, newname) {
		.Object@means <- renameReference(.Object@means, oldname, newname)
		.Object@covariance <- renameReference(.Object@covariance, oldname, newname)
		.Object@data <- renameReference(.Object@data, oldname, newname)
		.Object@thresholds <- sapply(.Object@thresholds, renameReference, oldname, newname)		
		return(.Object)
})

verifyExpectedObservedNames <- function(data, covName, flatModel, modelname, objectiveName) {
	covariance <- flatModel[[covName]]
	if (is(covariance, "MxMatrix") && !identical(dim(covariance), dim(data))) {
		msg <- paste("The dimensions for the expected covariance matrix",
			"and the observed covariance matrix",
			"in the", objectiveName, "expectation function in model",
			omxQuotes(modelname), "are not identical.")
		stop(msg, call. = FALSE)
	}
	if (!identical(dimnames(covariance), dimnames(data))) {
		msg <- paste("The dimnames for the expected covariance matrix",
			"and the observed covariance matrix",
			"in the", objectiveName, "expectation function in model",
			omxQuotes(modelname), "are not identical.")
		stop(msg, call. = FALSE)		
	}
}

verifyMeans <- function(meansName, mxDataObject, flatModel, modelname) {
	means <- flatModel[[meansName]]
	if (!is.null(means)) {
		if(any(is.na(mxDataObject@means))) {
			msg <- paste("In model", omxQuotes(modelname),
				"the Normal expectation function contains an expected means",
				"vector but the model is missing some data",
				"for the observed means.")
			stop(msg, call. = FALSE)
		}
		meanDimnames <- dimnames(means)		
	}
}

setMethod("genericExpFunConvert", "MxExpectationNormal", 
	function(.Object, flatModel, model, labelsData, defVars, dependencies) {
		modelname <- imxReverseIdentifier(model, .Object@name)[[1]]
		name <- .Object@name
		if(is.na(.Object@data)) {
			msg <- paste("The normal expectation function",
				"does not have a dataset associated with it in model",
				omxQuotes(modelname))
			stop(msg, call.=FALSE)
		}

		mxDataObject <- flatModel@datasets[[.Object@data]]
		dataName <- .Object@data
		.Object@data <- imxLocateIndex(flatModel, .Object@data, name)
		threshName <- .Object@thresholds
		.Object@thresholds <- imxLocateIndex(flatModel, threshName, name)
		covName <- .Object@covariance
		covariance <- flatModel[[covName]]
		.Object@covariance <- imxLocateIndex(flatModel, .Object@covariance, name)
		meansName <- .Object@means
		.Object@means <- imxLocateIndex(flatModel, .Object@means, name)

		if (inherits(mxDataObject, "MxDataDynamic")) return(.Object)

		if (mxDataObject@type != "raw") {
			verifyExpectedObservedNames(mxDataObject@observed, covName, flatModel, modelname, "Normal")
			verifyMeans(meansName, mxDataObject, flatModel, modelname)
		}
		verifyObservedNames(mxDataObject@observed, mxDataObject@means, mxDataObject@type, flatModel, modelname, "Normal")
		checkNumericData(mxDataObject)
		checkNumberOrdinalColumns(mxDataObject)
		covNames <- colnames(covariance)
		.Object@definitionVars <- imxFilterDefinitionVariables(defVars, dataName)
		verifyMvnNames(covName, meansName, "expected", flatModel, modelname, class(.Object))
		.Object@dataColumns <- generateDataColumns(flatModel, covNames, dataName)
		verifyThresholds(flatModel, model, labelsData, dataName, covNames, threshName)
		retval <- generateThresholdColumns(flatModel, model, labelsData, covNames, dataName, threshName)
		.Object@thresholdColumns <- retval[[1]] 
		.Object@thresholdLevels <- retval[[2]]
		if (single.na(.Object@dims)) {
			.Object@dims <- covNames
		}
		return(.Object)
})

verifyMvnNames <- function(covName, meansName, type, flatModel, modelname, expectationName) {
	if (is.na(meansName)) {
		means <- NA
	} else {
		means <- flatModel[[meansName]]
	}
	covariance <- flatModel[[covName]]
	if (length(covariance)) {
		covDimnames <- dimnames(covariance)
		if (is.null(covDimnames)) {
			msg <- paste("The",type,"covariance matrix associated",
				     "with", expectationName, "in model",
				     omxQuotes(modelname), "does not contain dimnames.")
			stop(msg, call. = FALSE)	
		}
		covRows <- covDimnames[[1]]
		covCols <- covDimnames[[2]]
		if (is.null(covRows) || is.null(covCols) ||
		    (length(covRows) != length(covCols)) || !all(covRows == covCols)) {
			msg <- paste("The",type,"covariance matrix associated",
				     "with", expectationName, "in model",
				     omxQuotes(modelname), "does not contain identical",
				     "row and column dimnames.")
			stop(msg, call.=FALSE)
		}
	}
	if (is.null(means) || (!isS4(means) && is.na(means)) || !length(means)) return()
	meanDimnames <- dimnames(means)
	if (is.null(meanDimnames)) {
			msg <- paste("The",type,"means matrix associated",
				"with", expectationName, "in model",
				omxQuotes(modelname), "does not contain dimnames.")
			stop(msg, call.=FALSE)	
	}
	meanRows <- meanDimnames[[1]]
	meanCols <- meanDimnames[[2]]
	meanNames <- c()
	if (length(meanRows) > length(meanCols)) { meanNames <- meanRows } else { meanNames <- meanCols }
	if ((length(covCols) != length(meanNames)) || !all(covCols == meanNames)) {
			msg <- paste("The",type,"covariance and",type,
				"means matrices associated",
				"with", expectationName, "in model",
				omxQuotes(modelname), "do not contain identical",
				"dimnames.")
			stop(msg, call.=FALSE)
	}
}

verifyObservedNames <- function(data, means, type, flatModel, modelname, expectationName) {
	dataNames <- dimnames(data)
	if(is.null(dataNames)) {
		msg <- paste("The observed data associated with the",
			expectationName, "expectation function in model",
			omxQuotes(modelname), "does not contain dimnames.")
		stop(msg, call. = FALSE)
	}
	if (type == "cov" || type == "cor") {
		if (length(dataNames) < 2 ||
			is.null(dataNames[[1]]) || is.null(dataNames[[2]]) || 
			!identical(dataNames[[1]], dataNames[[2]])) {
				msg <- paste("The dataset associated with the", expectationName,
					"expectation function in model", omxQuotes(modelname),
    	            "does not contain identical row and column non-NULL dimnames.")
			stop(msg, call. = FALSE)
		}
		if (!single.na(means) && is.null(dimnames(means))) {
			msg <- paste("In model", omxQuotes(modelname), 
				", the observed means vector does not contain column names.",
				"Use the name() function to assign names to the means vector.")
			stop(msg, call. = FALSE)
		}
		if (!single.na(means) && !identical(dataNames[[1]], dimnames(means)[[2]])) {
			msg <- paste("The observed covariance or correlation matrix associated with the", expectationName,
				"expectation function in model", omxQuotes(modelname),
				"does not contain identical dimnames to the observed means vector.")
			stop(msg, call. = FALSE)
		}
	} else if ((type == "raw") && (length(dataNames) < 2 || is.null(dataNames[[2]]))) {
		msg <- paste("The dataset associated with the", expectationName,
				"expectation function in model", omxQuotes(modelname),
				"does not contain column names (use dimnames).")
		stop(msg, call. = FALSE)
	}
}

generateDataColumns <- function(flatModel, covNames, dataName) {
	retval <- c()
	dataColumnNames <- dimnames(flatModel@datasets[[dataName]]@observed)[[2]]
	for(i in 1:length(covNames)) {
		targetName <- covNames[[i]]
		index <- match(targetName, dataColumnNames)
		if(is.na(index)) {
			msg <- paste("The column name", omxQuotes(targetName),
				"in the expected covariance matrix",
				"of the expectation function in model",
				omxQuotes(flatModel@name),
				"cannot be found in the column names of the data.")
			stop(msg, call. = FALSE)
		}
		retval[[i]] <- index - 1
	}
	return(retval)
}


updateThresholdDimnames <- function(flatExpectation, flatModel, labelsData) {
	threshName <- flatExpectation@thresholds
	if (is.na(threshName)) {
		return(flatModel)
	}
	thresholds <- flatModel[[threshName]]
	if (is.null(thresholds)) {
		modelname <- getModelName(flatExpectation)
		stop(paste("Unknown thresholds name", 
			omxQuotes(simplifyName(threshName, modelname)),
			"detected in the expectation function",
			"of model", omxQuotes(modelname)), call. = FALSE)
	}
	dims <- flatExpectation@threshnames
	if (!is.null(colnames(thresholds)) && !single.na(dims) && 
		!identical(colnames(thresholds), dims)) {
		modelname <- getModelName(flatExpectation)
		msg <- paste("The thresholds matrix associated",
		"with the expectation function in model", 
		omxQuotes(modelname), "contains column names and",
		"the expectation function has specified non-identical threshnames.")
		stop(msg, call.=FALSE)      
	}
	if (is.null(colnames(thresholds)) && !single.na(dims)) {
		tuple <- evaluateMxObject(threshName, flatModel, labelsData, new.env(parent = emptyenv()))
		threshMatrix <- tuple[[1]]
		if (ncol(threshMatrix) != length(dims)) {
			modelname <- getModelName(flatExpectation)
			msg <- paste("The thresholds matrix associated",
			"with the expectation function in model", 
			omxQuotes(modelname), "is not of the same length as the 'threshnames'",
			"argument provided by the expectation function. The 'threshnames' argument is",
			"of length", length(dims), "and the expected covariance matrix",
			"has", ncol(threshMatrix), "columns.")
			stop(msg, call.=FALSE)      
		}
		dimnames(flatModel[[threshName]]) <- list(NULL, dims)
	}
	return(flatModel)
}

updateExpectationDimnames <- function(flatExpectation, flatModel,
		labelsData, unsafe = FALSE) {
	covName <- flatExpectation@covariance
	meansName <- flatExpectation@means
	if (is.na(meansName)) {
		means <- NA
	} else {
		means <- flatModel[[meansName]]
	}
	covariance <- flatModel[[covName]]
	if (is.null(covariance)) {
		modelname <- getModelName(flatExpectation)
		stop(paste("Unknown expected covariance name", 
			omxQuotes(simplifyName(covName, modelname)),
			"detected in the expectation function",
			"of model", omxQuotes(modelname)), call. = FALSE)
	}
	if (is.null(means)) {
		modelname <- getModelName(flatExpectation)
		stop(paste("Unknown expected means name", 
			omxQuotes(simplifyName(meansName, modelname)),
			"detected in the expectation function",
			"of model", omxQuotes(modelname)), call. = FALSE)
	}
	dims <- flatExpectation@dims
	if (!is.null(dimnames(covariance)) && !single.na(dims) && 
		!identical(dimnames(covariance), list(dims, dims))) {
		modelname <- getModelName(flatExpectation)
		msg <- paste("The expected covariance matrix associated",
			"with the expectation function in model", 
			omxQuotes(modelname), "contains dimnames: ", 
            paste(toString(dimnames(covariance)), ".", sep = ""),
			"The expectation function has specified dimnames:", 
			paste(toString(dims), ".", sep =""))
		stop(msg, call.=FALSE)		
	}
	if (is.null(dimnames(covariance)) && !single.na(dims)) {
		if (!unsafe) {
			tuple <- evaluateMxObject(covName, flatModel, labelsData, new.env(parent = emptyenv()))
			covMatrix <- tuple[[1]]
			if (nrow(covMatrix) != ncol(covMatrix)) {
				modelname <- getModelName(flatExpectation)
				msg <- paste("The expected covariance matrix associated",
					"with the expectation function in model", 
					omxQuotes(modelname), "is not a square matrix.")
				stop(msg, call.=FALSE)		
			}
			if (nrow(covMatrix) != length(dims)) {
				modelname <- getModelName(flatExpectation)
				msg <- paste("The expected covariance matrix associated",
					"with the expectation function in model", 
					omxQuotes(modelname), "is not of the same length as the 'dimnames'",
					"argument provided by the expectation function. The 'dimnames' argument is",
					"of length", length(dims), "and the expected covariance matrix",
					"has", nrow(covMatrix), "rows and columns.")
				stop(msg, call.=FALSE)		
			}
		}
		dimnames(flatModel[[covName]]) <- list(dims, dims)
	}

	if (!isS4(means) && is.na(means)) {
		return(flatModel)
	}

	if (!is.null(dimnames(means)) && !single.na(dims) &&
		!identical(dimnames(means), list(NULL, dims))) {
		modelname <- getModelName(flatExpectation)
		msg <- paste("The expected means matrix associated",
			"with the expectation function in model", 
			omxQuotes(modelname), "contains dimnames: ", 
            paste(toString(dimnames(means)), ".", sep = ""),
			"The expectation function has specified dimnames:", 
			paste(toString(dims), ".", sep =""))
		stop(msg, call.=FALSE)	
	}
	if (is.null(dimnames(means)) && !single.na(dims)) {
		if (!unsafe) {
			tuple <- evaluateMxObject(meansName, flatModel, labelsData, new.env(parent = emptyenv()))
			meansMatrix <- tuple[[1]]
			if (nrow(meansMatrix) != 1) {
				modelname <- getModelName(flatExpectation)
				msg <- paste("The expected means vector associated",
					"with the expectation function in model", 
					omxQuotes(modelname), "is not a 1 x n matrix.",
					"It has dimensions", nrow(meansMatrix), "x", 
					paste(ncol(meansMatrix), '.', sep=''))
				stop(msg, call.=FALSE)		
			}
			if (ncol(meansMatrix) != length(dims)) {
				modelname <- getModelName(flatExpectation)
				msg <- paste("The expected means vector associated",
					"with the expectation function in model", 
					omxQuotes(modelname), "is not of the same length as the 'dimnames'",
					"argument provided by the expectation function. The 'dimnames' argument is",
					"of length", length(dims), "and the expected means vector",
					"has", ncol(meansMatrix), "columns.")
				stop(msg, call.=FALSE)
			}
		}
		dimnames(flatModel[[meansName]]) <- list(NULL, dims)
	}
	return(flatModel)
}



mxExpectationNormal <- function(covariance, means = NA, 
	dimnames = NA, thresholds = NA, threshnames = dimnames) {
	if (missing(covariance) || typeof(covariance) != "character") {
		stop("'covariance' argument is not a string (the name of the expected covariance matrix)")
	}
	if (!(single.na(means) || typeof(means) == "character")) {
		stop("Means argument is not a string (the name of the expected means matrix)")
	}
	if (is.na(means)) means <- as.integer(NA)
	if (single.na(thresholds)) thresholds <- as.character(NA)
	if (single.na(dimnames)) dimnames <- as.character(NA)
	if (single.na(threshnames)) threshnames <- as.character(NA)
	if (!is.vector(dimnames) || typeof(dimnames) != 'character') {
		stop("'dimnames' argument is not a character vector")
	}
	if (!is.vector(threshnames) || typeof(threshnames) != 'character') {
		stop("'threshnames' argument is not a character vector")
	}
	if (length(thresholds) != 1) {
		stop("'thresholds' argument must be a single matrix or algebra name")
	}
	if (length(dimnames) == 0) {
		stop("'dimnames' argument cannot be an empty vector")
	}
	if (length(threshnames) == 0) {
		stop("'threshnames' argument cannot be an empty vector")
	}
	if (length(dimnames) > 1 && any(is.na(dimnames))) {
		stop("NA values are not allowed for 'dimnames' vector")
	}
	if (length(threshnames) > 1 && any(is.na(threshnames))) {
		stop("NA values are not allowed for 'threshnames' vector")
	}
	return(new("MxExpectationNormal", covariance, means, dimnames, thresholds, threshnames))
}

displayMxExpectationNormal <- function(expectation) {
	cat("MxExpectationNormal", omxQuotes(expectation@name), '\n')
	cat("$covariance :", omxQuotes(expectation@covariance), '\n')
	cat("$means :", omxQuotes(expectation@means), '\n')
	if (single.na(expectation@dims)) {
		cat("$dims : NA \n")
	} else {
		cat("$dims :", omxQuotes(expectation@dims), '\n')
	}
	if (single.na(expectation@thresholds)) {
		cat("$thresholds : NA \n")
	} else {
		cat("$thresholds :", omxQuotes(expectation@thresholds), '\n')
	}
	if (single.na(expectation@threshnames)) {
		cat("$threshnames : NA \n")
	} else {
		cat("$threshnames :", omxQuotes(expectation@threshnames), '\n')
	}
	invisible(expectation)
}


setMethod("print", "MxExpectationNormal", function(x,...) { 
	displayMxExpectationNormal(x) 
})

setMethod("show", "MxExpectationNormal", function(object) { 
	displayMxExpectationNormal(object) 
})
