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


setClass(Class = "MxFIMLObjective",
	representation = representation(
		covariance = "MxCharOrNumber",
		means = "MxCharOrNumber",
		definitionVars = "list",
		thresholds = "MxCharOrNumber",
		dims = "character",
		dataColumns = "numeric",
		thresholdColumns = "list",
		vector = "logical",
		threshnames = "character",
		metadata = "MxBaseObjectiveMetaData"),
	contains = "MxBaseObjective")

setMethod("initialize", "MxFIMLObjective",
	function(.Object, covariance, means, dims, thresholds, vector, threshnames,
		data = as.integer(NA), definitionVars = list(), name = 'objective') {
		.Object@name <- name
		.Object@covariance <- covariance
		.Object@means <- means
		.Object@data <- data
		.Object@definitionVars <- definitionVars
		.Object@thresholds <- thresholds
		.Object@dims <- dims
		.Object@vector <- vector
		.Object@threshnames <- threshnames
		return(.Object)
	}
)

setMethod("genericObjDependencies", signature("MxFIMLObjective"),
	function(.Object, dependencies) {
	sources <- c(.Object@covariance, .Object@means, .Object@thresholds)
	sources <- sources[!is.na(sources)]
	dependencies <- imxAddDependency(sources, .Object@name, dependencies)
	return(dependencies)
})


setMethod("genericObjFunNamespace", signature("MxFIMLObjective"), 
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
		if (!is.null(.Object@metadata)) {
			metadata <- .Object@metadata
			metadata@A <-  imxConvertIdentifier(metadata@A, modelname, namespace)
			metadata@S <-  imxConvertIdentifier(metadata@S, modelname, namespace)
			metadata@F <-  imxConvertIdentifier(metadata@F, modelname, namespace)
			metadata@M <-  imxConvertIdentifier(metadata@M, modelname, namespace)
			.Object@metadata <- metadata
		}
		return(.Object)
})

setMethod("genericObjRename", signature("MxFIMLObjective"),
	function(.Object, oldname, newname) {
		.Object@covariance <- renameReference(.Object@covariance, oldname, newname)
		.Object@means <- renameReference(.Object@means, oldname, newname)		
		.Object@data <- renameReference(.Object@data, oldname, newname)	
		.Object@thresholds <- sapply(.Object@thresholds, renameReference, oldname, newname)		
		return(.Object)
})

setMethod("genericObjFunConvert", signature("MxFIMLObjective"), 
	function(.Object, flatModel, model, defVars) {
		modelname <- imxReverseIdentifier(model, .Object@name)[[1]]
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
			stop(msg, call. = FALSE)
		}
		verifyObservedNames(mxDataObject@observed, mxDataObject@type, flatModel, modelname, "FIML")
		checkNumericData(mxDataObject)
		meansName <- .Object@means
		covName <- .Object@covariance
		dataName <- .Object@data
		threshName <- .Object@thresholds
		.Object@definitionVars <- imxFilterDefinitionVariables(defVars, .Object@data)
		.Object@means <- imxLocateIndex(flatModel, .Object@means, name)
		.Object@covariance <- imxLocateIndex(flatModel, .Object@covariance, name)
		.Object@data <- imxLocateIndex(flatModel, .Object@data, name)
		verifyExpectedNames(covName, meansName, flatModel, modelname, "FIML")
		.Object@dataColumns <- generateDataColumns(flatModel, covName, dataName)
		verifyThresholds(flatModel, model, dataName, covName, threshName)
		.Object@thresholds <- imxLocateIndex(flatModel, threshName, name)
		if (!is.null(.Object@metadata)) {
			metadata <- .Object@metadata
			metadata@A <-  imxLocateIndex(flatModel, metadata@A, name)
			metadata@S <-  imxLocateIndex(flatModel, metadata@S, name)
			metadata@F <-  imxLocateIndex(flatModel, metadata@F, name)
			metadata@M <-  imxLocateIndex(flatModel, metadata@M, name)
			.Object@metadata <- metadata
		}
		.Object@thresholdColumns <- generateThresholdColumns(flatModel, 
			model, dataName, threshName)
		if (length(mxDataObject@observed) == 0) {
			.Object@data <- as.integer(NA)
		}
		if (single.na(.Object@dims)) {
			.Object@dims <- rownames(flatModel[[covName]])
		}
		return(.Object)
})

setMethod("genericObjModelConvert", "MxFIMLObjective",
	function(.Object, job, model, namespace, flatJob) {
		job <- updateObjectiveDimnames(.Object, job, model@name, "FIML")
		job <- updateThresholdDimnames(.Object, job, model@name)
		precision <- "Function Precision"
		if(!single.na(.Object@thresholds) && 
			is.null(job@options[[precision]])) {
			job <- mxOption(job, precision, 1e-9)
		}
		job@.newobjects <- TRUE
		job@.newobjective <- FALSE
		job@.newtree <- FALSE
		return(job)
	}
)

setMethod("genericObjInitialMatrix", "MxFIMLObjective",
	function(.Object, flatModel) {
		flatObjective <- flatModel@objectives[[.Object@name]]
		if (flatObjective@vector == FALSE) {
			return(matrix(as.double(NA), 1, 1))
		} else {
			modelname <- imxReverseIdentifier(flatModel, flatObjective@name)[[1]]
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

updateThresholdDimnames <- function(flatObjective, job, modelname) {
	threshName <- flatObjective@thresholds
	if (is.na(threshName)) {
		return(job)
	}
	thresholds <- job[[threshName]]
	if (is.null(thresholds)) {
		stop(paste("Unknown thresholds name", 
			omxQuotes(simplifyName(thresholds, modelname)),
			"detected in the objective function",
			"of model", omxQuotes(modelname)), call. = FALSE)
	}
	dims <- flatObjective@threshnames
	if (!is.null(colnames(thresholds)) && !single.na(dims) && 
		!identical(colnames(thresholds), dims)) {
		msg <- paste("The thresholds matrix associated",
			"with the FIML objective in model", 
			omxQuotes(modelname), "contains column names and",
			"the objective function has specified non-identical threshnames.")
		stop(msg, call.=FALSE)		
	}
	if (is.null(colnames(thresholds)) && !single.na(dims)) {
		threshMatrix <- eval(substitute(mxEval(x, job, compute=TRUE), list(x = as.symbol(threshName))))
		if (ncol(threshMatrix) != length(dims)) {
			msg <- paste("The thresholds matrix associated",
				"with the FIML objective in model", 
				omxQuotes(modelname), "is not of the same length as the 'threshnames'",
				"argument provided by the objective function. The 'threshnames' argument is",
				"of length", length(dims), "and the expected covariance matrix",
				"has", ncol(threshMatrix), "columns.")
			stop(msg, call.=FALSE)		
		}
		dimnames(thresholds) <- list(NULL, dims)
		job[[threshName]] <- thresholds
	}
	return(job)
}

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
			omxQuotes(simplifyName(covName, modelname)),
			"detected in the objective function",
			"of model", omxQuotes(modelname)), call. = FALSE)
	}
	if (is.null(means)) {
		stop(paste("Unknown expected means name", 
			omxQuotes(simplifyName(meansName, modelname)),
			"detected in the objective function",
			"of model", omxQuotes(modelname)), call. = FALSE)
	}
	dims <- flatObjective@dims
	if (!is.null(dimnames(covariance)) && !single.na(dims) && 
		!identical(dimnames(covariance), list(dims, dims))) {
		msg <- paste("The expected covariance matrix associated",
			"with the", objectiveName, "objective in model", 
			omxQuotes(modelname), "contains dimnames and",
			"the objective function has specified dimnames.")
		stop(msg, call.=FALSE)		
	}
	if (is.null(dimnames(covariance)) && !single.na(dims)) {
		covMatrix <- eval(substitute(mxEval(x, job, compute=TRUE), list(x = as.symbol(covName))))
		if (nrow(covMatrix) != ncol(covMatrix)) {
			msg <- paste("The expected covariance matrix associated",
				"with the", objectiveName, "objective in model", 
				omxQuotes(modelname), "is not a square matrix.")
			stop(msg, call.=FALSE)		
		}
		if (nrow(covMatrix) != length(dims)) {
			msg <- paste("The expected covariance matrix associated",
				"with the", objectiveName, "objective in model", 
				omxQuotes(modelname), "is not of the same length as the 'dimnames'",
				"argument provided by the objective function. The 'dimnames' argument is",
				"of length", length(dims), "and the expected covariance matrix",
				"has", nrow(covMatrix), "rows and columns.")
			stop(msg, call.=FALSE)		
		}
		dimnames(covariance) <- list(dims, dims)
		job[[covName]] <- covariance
	}
	if (!isS4(means) && is.na(means)) return(job)
	if (!is.null(dimnames(means)) && !single.na(dims) &&
		!identical(dimnames(means), list(NULL, dims))) {
		msg <- paste("The expected means matrix associated",
			"with the", objectiveName, "objective in model", 
			omxQuotes(modelname), "contains dimnames and",
			"the objective function has specified dimnames.")
		stop(msg, call.=FALSE)	
	}
	if (is.null(dimnames(means)) && !single.na(dims)) {
		meansMatrix <- eval(substitute(mxEval(x, job, compute=TRUE), 
			list(x = as.symbol(meansName))))
		if (nrow(meansMatrix) != 1) {
			msg <- paste("The expected means vector associated",
				"with the", objectiveName, "objective in model", 
				omxQuotes(modelname), "is not a 1 x n matrix.",
				"It has dimensions", nrow(meansMatrix), "x", 
				paste(ncol(meansMatrix), '.', sep=''))
			stop(msg, call.=FALSE)		
		}
		if (ncol(meansMatrix) != length(dims)) {
			msg <- paste("The expected means vector associated",
				"with the", objectiveName, "objective in model", 
				omxQuotes(modelname), "is not of the same length as the 'dimnames'",
				"argument provided by the objective function. The 'dimnames' argument is",
				"of length", length(dims), "and the expected means vector",
				"has", ncol(meansMatrix), "columns.")
			stop(msg, call.=FALSE)
		}
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
	meanDimnames <- dimnames(means)
	if (is.null(meanDimnames)) {
			msg <- paste("The expected means matrix associated",
				"with the", objectiveName, "objective function in model", 
				omxQuotes(modelname), "does not contain dimnames.")
			stop(msg, call.=FALSE)	
	}
	meanRows <- meanDimnames[[1]]
	meanCols <- meanDimnames[[2]]
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
				"in the expected covariance matrix",
				"of the FIML objective function in model",
				omxQuotes(flatModel@name),
				"cannot be found in the column names of the data.")
			stop(msg, call. = FALSE)
		}
		retval[[i]] <- index - 1
	}
	return(retval)
}


mxFIMLObjective <- function(covariance, means, dimnames = NA, 
						thresholds = NA, vector = FALSE, threshnames = dimnames) {
	if (missing(covariance) || typeof(covariance) != "character") {
		stop("'covariance' argument is not a string (the name of the expected covariance matrix)")
	}
	if (missing(means) || typeof(means) != "character") {
		stop("'means' argument is not a string (the name of the expected means vector)")
	}
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
	if (length(vector) > 1 || typeof(vector) != "logical") {
		stop("'vector' argument is not a logical value")
	}
	return(new("MxFIMLObjective", covariance, means, dimnames, thresholds, vector, threshnames))
}

displayFIMLObjective <- function(objective) {
	cat("MxFIMLObjective", omxQuotes(objective@name), '\n')
	cat("@covariance :", omxQuotes(objective@covariance), '\n')
	cat("@means :", omxQuotes(objective@means), '\n')
	cat("@vector :", objective@vector, '\n')
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
	if (single.na(objective@threshnames)) {
		cat("@threshnames : NA \n")
	} else {
		cat("@threshnames :", omxQuotes(objective@threshnames), '\n')
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
