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


setClass(Class = "MxFIMLObjective",
	representation = representation(
		covariance = "MxCharOrNumber",
		means = "MxCharOrNumber",
		definitionVars = "list",
		thresholds = "character",
		dims = "character",
		dataColumns = "numeric",
		thresholdColumns = "list",
		vector = "logical"),
	contains = "MxBaseObjective")

setMethod("initialize", "MxFIMLObjective",
	function(.Object, covariance, means, dims, thresholds, vector, 
		data = as.integer(NA), definitionVars = list(), name = 'objective') {
		.Object@name <- name
		.Object@covariance <- covariance
		.Object@means <- means
		.Object@data <- data
		.Object@definitionVars <- definitionVars
		.Object@thresholds <- thresholds
		.Object@dims <- dims
		.Object@vector <- vector
		return(.Object)
	}
)

setMethod("omxObjDependencies", signature("MxFIMLObjective"),
	function(.Object, dependencies) {
	sources <- c(.Object@covariance, .Object@means, .Object@thresholds)
	sources <- sources[!is.na(sources)]
	dependencies <- omxAddDependency(sources, .Object@name, dependencies)
	return(dependencies)
})


setMethod("omxObjFunNamespace", signature("MxFIMLObjective"), 
	function(.Object, modelname, namespace) {
		.Object@name <- omxIdentifier(modelname, .Object@name)
		.Object@covariance <- omxConvertIdentifier(.Object@covariance, 
			modelname, namespace)
		.Object@means <- omxConvertIdentifier(.Object@means, 
			modelname, namespace)
		.Object@data <- omxConvertIdentifier(.Object@data, 
			modelname, namespace)
		.Object@thresholds <- sapply(.Object@thresholds, function(x) {omxConvertIdentifier(x,
			modelname, namespace)})
		return(.Object)
})

setMethod("omxObjFunConvert", signature("MxFIMLObjective"), 
	function(.Object, flatModel, model) {
		modelname <- omxReverseIdentifier(model, .Object@name)[[1]]
		name <- .Object@name
		if(is.na(.Object@data)) {
			msg <- paste("The FIML objective",
				"does not have a dataset associated with it in model",
				omxQuotes(modelname))
			stop(msg, call.=FALSE)
		}
		mxDataObject <- flatModel@datasets[[.Object@data]]
		if (mxDataObject@type != 'raw') {
			msg <- paste("The dataset associated with the FIML objective", 
				"in model", omxQuotes(modelname), "is not raw data.")
			stop(msg, call.=FALSE)
		}
		verifyObservedNames(mxDataObject@observed, mxDataObject@type, flatModel, modelname, "FIML")
		checkNumericData(mxDataObject)
		meansName <- .Object@means
		covName <- .Object@covariance
		dataName <- .Object@data
		threshNames <- .Object@thresholds
		.Object@means <- omxLocateIndex(flatModel, .Object@means, name)
		.Object@covariance <- omxLocateIndex(flatModel, .Object@covariance, name)
		.Object@data <- omxLocateIndex(flatModel, .Object@data, name)
		verifyExpectedNames(covName, meansName, flatModel, modelname, "FIML")
		.Object@definitionVars <- generateDefinitionList(flatModel)
		.Object@dataColumns <- generateDataColumns(flatModel, covName, dataName)
		.Object@thresholdColumns <- generateThresholdColumns(flatModel, threshNames, dataName)
		if (length(mxDataObject@observed) == 0) {
			.Object@data <- as.integer(NA)
		}
		return(.Object)
})

setMethod("omxObjModelConvert", "MxFIMLObjective",
	function(.Object, job, model, flatJob) {
		job <- updateObjectiveDimnames(.Object, job, model@name, "FIML")
		return(job)
	}
)

setMethod("omxObjInitialMatrix", "MxFIMLObjective",
	function(.Object, flatModel) {
		flatObjective <- flatModel@objectives[[.Object@name]]
		if (flatObjective@vector == FALSE) {
			return(matrix(as.double(NA), 1, 1))
		} else {
			modelname <- omxReverseIdentifier(flatModel, flatObjective@name)[[1]]
			name <- flatObjective@name
			if(is.na(flatObjective@data)) {
				msg <- paste("The FIML objective",
				"does not have a dataset associated with it in model",
				omxQuotes(modelname))
				stop(msg, call.=FALSE)
			}
			mxDataObject <- flatModel@datasets[[flatObjective@data]]
			if (mxDataObject@type != 'raw') {
				msg <- paste("The dataset associated with the FIML objective", 
					"in model", omxQuotes(modelname), "is not raw data.")
				stop(msg, call.=FALSE)
			}
			rows <- nrow(mxDataObject@observed)
			return(matrix(as.double(NA), rows, 1))
		}
})


updateObjectiveDimnames <- function(flatObjective, job, modelname, objectiveName) {
	covName <- flatObjective@covariance
	meansName <- flatObjective@means
	if (is.na(meansName)) {
		means <- NA
	} else {
		means <- job[[meansName]]
	}
	covariance <- job[[covName]]
	if (is.null(covariance)) {
		stop(paste("Unknown expected covariance name", 
			omxQuotes(simplifyName(covName)),
			"detected in the objective function",
			"of model", omxQuotes(modelname)), call. = FALSE)
	}
	if (is.null(means)) {
		stop(paste("Unknown expected means name", 
			omxQuotes(simplifyName(meansName)),
			"detected in the objective function",
			"of model", omxQuotes(modelname)), call. = FALSE)
	}
	dims <- flatObjective@dims
	if (!is.null(dimnames(covariance)) && !single.na(dims) && 
		!identical(dimnames(covariance), list(dims, dims))) {
		msg <- paste("The expected covariance matrix associated",
			"with the", objectiveName, "objective in model", 
			omxQuotes(modelname), "contains dimnames and",
			"the objective function has specified dimnames")
		stop(msg, call.=FALSE)		
	}
	if (is.null(dimnames(covariance)) && !single.na(dims)) {
		dimnames(covariance) <- list(dims, dims)
		job[[covName]] <- covariance
	}
	if (!isS4(means) && is.na(means)) return(job)
	if (!is.null(dimnames(means)) && !single.na(dims) &&
		!identical(dimnames(means), list(NULL, dims))) {
		msg <- paste("The expected means matrix associated",
			"with the", objectiveName, "objective in model", 
			omxQuotes(modelname), "contains dimnames and",
			"the objective function has specified dimnames")
		stop(msg, call.=FALSE)	
	}
	if (is.null(dimnames(means)) && !single.na(dims)) {
		dimnames(means) <- list(NULL, dims)
		job[[meansName]] <- means
	}
	return(job)
}

verifyObservedNames <- function(data, type, flatModel, modelname, objectiveName) {
	dataNames <- dimnames(data)
	if(is.null(dataNames)) {
		msg <- paste("The observed data associated with the",
			objectiveName, "objective in model",
			omxQuotes(modelname), "does not contain dimnames.")
		stop(msg, call. = FALSE)
	}
	if ((type == "cov" || type == "cor") && (length(dataNames) < 2 ||
		is.null(dataNames[[1]]) || is.null(dataNames[[2]]) || 
		!identical(dataNames[[1]], dataNames[[2]]))) {
                msg <- paste("The dataset associated with the", objectiveName,
                                "objective in model", omxQuotes(modelname),
                                "does not contain identical row and column non-NULL dimnames.")
                stop(msg, call. = FALSE)
	} else if ((type == "raw") && (length(dataNames) < 2 || is.null(dataNames[[2]]))) {
		msg <- paste("The dataset associated with the", objectiveName,
				"objective in model", omxQuotes(modelname),
				"does not contain column names (use dimnames).")
		stop(msg, call. = FALSE)
	}
}

verifyExpectedNames <- function(covName, meansName, flatModel, modelname, objectiveName) {
	if (is.na(meansName)) {
		means <- NA
	} else {
		means <- flatModel[[meansName]]
	}
	covariance <- flatModel[[covName]]
	covariance <- dimnames(covariance)
	if (is.null(covariance)) {
			msg <- paste("The expected covariance matrix associated",
				"with the", objectiveName, "objective in model", 
				omxQuotes(modelname), "does not contain dimnames.")
			stop(msg, call. = FALSE)	
	}
	covRows <- covariance[[1]]
	covCols <- covariance[[2]]	
	if (is.null(covRows) || is.null(covCols) ||
		(length(covRows) != length(covCols)) || !all(covRows == covCols)) {
			msg <- paste("The expected covariance matrix associated",
				"with the", objectiveName, "objective in model", 
				omxQuotes(modelname), "does not contain identical",
				"row and column dimnames.")
			stop(msg, call.=FALSE)
	}
	if (!isS4(means) && is.na(means)) return()
	means <- dimnames(means)
	if (is.null(means)) {
			msg <- paste("The expected means matrix associated",
				"with the", objectiveName, "objective function in model", 
				omxQuotes(modelname), "does not contain dimnames.")
			stop(msg, call.=FALSE)	
	}
	meanRows <- means[[1]]
	meanCols <- means[[2]]
	if (!is.null(meanRows) && length(meanRows) > 1) {
		msg <- paste("The expected means matrix associated",
			"with the", objectiveName, "objective in model", 
			omxQuotes(modelname), "is not a 1 x N matrix.")
			stop(msg, call.=FALSE)
	}
	if ((length(covCols) != length(meanCols)) || !all(covCols == meanCols)) {
			msg <- paste("The expected covariance and expected",
				"means matrices associated",
				"with the", objectiveName, "objective function in model", 
				omxQuotes(modelname), "do not contain identical",
				"dimnames.")
			stop(msg, call.=FALSE)
	}
}

generateDataColumns <- function(flatModel, covName, dataName) {
	retval <- c()
	definitionNames <- dimnames(flatModel@datasets[[dataName]]@observed)[[2]]
	covariance <- flatModel[[covName]]
	covNames <- dimnames(covariance)[[2]]
	for(i in 1:length(covNames)) {
		targetName <- covNames[[i]]
		index <- match(targetName, definitionNames)
		if(is.na(index)) {
			msg <- paste("The column name", omxQuotes(targetName),
				"in the observed covariance matrix",
				"of the FIML objective function in model",
				omxQuotes(flatModel@name),
				"cannot be found in the dimnames of the data.")
			stop(msg, call. = FALSE)
		}
		retval[[i]] <- index - 1
	}
	return(retval)
}


mxFIMLObjective <- function(covariance, means, dimnames = NA, thresholds = NA, vector = FALSE) {
	if (missing(covariance) || typeof(covariance) != "character") {
		stop("Covariance argument is not a string (the name of the expected covariance matrix)")
	}
	if (missing(means) || typeof(means) != "character") {
		stop("Means argument is not a string (the name of the expected means vector)")
	}
	if (all(is.na(thresholds))) thresholds <- as.character(NA)
	if (single.na(dimnames)) dimnames <- as.character(NA)
	if (!is.vector(dimnames) || typeof(dimnames) != 'character') {
		stop("Dimnames argument is not a character vector")
	}
	if (length(dimnames) == 0) {
		stop("Dimnames argument cannot be an empty vector")
	}
	if (length(dimnames) > 1 && any(is.na(dimnames))) {
		stop("NA values are not allowed for dimnames vector")
	}
	if (length(vector) > 1 || typeof(vector) != "logical") {
		stop("Vector argument is not a logical value (return a vector or a scalar)")
	}	
	return(new("MxFIMLObjective", covariance, means, dimnames, thresholds, vector))
}

displayFIMLObjective <- function(objective) {
	cat("MxFIMLObjective", omxQuotes(objective@name), '\n')
	cat("@covariance :", omxQuotes(objective@covariance), '\n')
	cat("@means :", omxQuotes(objective@means), '\n')
	if (single.na(objective@dims)) {
		cat("@dims : NA \n")
	} else {
		cat("@dims :", omxQuotes(objective@dims), '\n')
	}	
	if (single.na(objective@thresholds)) {
		cat("@thresholds : NA \n")
	} else {
		cat("@thresholds :", omxQuotes(objective@thresholds), '\n')
	}	
	if (length(objective@result) == 0) {
		cat("@result: (not yet computed) ")
	} else {
		cat("@result:\n")
	}
	print(objective@result)
	invisible(objective)
}


setMethod("print", "MxFIMLObjective", function(x,...) { 
	displayFIMLObjective(x) 
})

setMethod("show", "MxFIMLObjective", function(object) { 
	displayFIMLObjective(object) 
})
