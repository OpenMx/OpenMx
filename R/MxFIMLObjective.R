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
		dataRow = "numeric"),
	contains = "MxBaseObjective")

setMethod("initialize", "MxFIMLObjective",
	function(.Object, covariance, means, data = as.numeric(NA),
		definitionVars = list(), name = 'objective') {
		.Object@name <- name
		.Object@covariance <- covariance
		.Object@means <- means
		.Object@data <- data
		.Object@definitionVars <- definitionVars
		return(.Object)
	}
)

setMethod("omxObjFunNamespace", signature("MxFIMLObjective"), 
	function(.Object, modelname, namespace) {
		.Object@name <- omxIdentifier(modelname, .Object@name)
		.Object@covariance <- omxConvertIdentifier(.Object@covariance, 
			modelname, namespace)
		.Object@means <- omxConvertIdentifier(.Object@means, 
			modelname, namespace)
		.Object@data <- omxConvertIdentifier(.Object@data, 
			modelname, namespace)
		return(.Object)
})

setMethod("omxObjFunConvert", signature("MxFIMLObjective"), 
	function(.Object, flatModel, model) {
		name <- .Object@name
		if(is.na(.Object@data)) {
			msg <- paste("The FIML objective",
				"does not have a dataset associated with it in model",
				omxQuotes(flatModel@name))
			stop(msg, call.=FALSE)
		}
		if (flatModel@datasets[[.Object@data]]@type != 'raw') {
			msg <- paste("The dataset associated with the FIML objective", 
				"in model", omxQuotes(flatModel@name), "is not raw data.")
			stop(msg, call.=FALSE)
		}
		dataNames <- dimnames(flatModel@datasets[[.Object@data]]@matrix)
		if (is.null(dataNames)) {
			msg <- paste("The dataset associated with the FIML objective", 
				"in model", omxQuotes(flatModel@name), 
				"does not contain column names (use dimnames).")
			stop(msg, call.=FALSE)
		}
		columnNames <- dataNames[[2]]
		meansName <- .Object@means
		covName <- .Object@covariance
		dataName <- .Object@data
		if(!is.na(.Object@means)) {			
			.Object@means <- omxLocateIndex(flatModel, .Object@means, name)
		}
		.Object@covariance <- omxLocateIndex(flatModel, .Object@covariance, name)
		.Object@data <- omxLocateIndex(flatModel, .Object@data, name)
		verifyExpectedNames(covName, meansName, flatModel)
		.Object@definitionVars <- generateDefinitionList(flatModel)
		.Object@dataRow <- generateDataRow(flatModel, covName, dataName)
		return(.Object)
})

verifyExpectedNames <- function(covName, meansName, flatModel) {
	if (is.na(meansName)) {
		means <- NA
	} else {
		means <- flatModel[[meansName]]
	}
	covariance <- flatModel[[covName]]
	covariance <- dimnames(covariance)
	if (is.null(covariance)) {
			msg <- paste("The expected covariance matrix associated",
				"with the FIML objective in model", 
				omxQuotes(flatModel@name), "does not contain dimnames.")
			stop(msg, call.=FALSE)	
	}
	covRows <- covariance[[1]]
	covCols <- covariance[[2]]	
	if (is.null(covRows) || is.null(covCols) ||
		(length(covRows) != length(covCols)) || !all(covRows == covCols)) {
			msg <- paste("The expected covariance matrix associated",
				"with the FIML objective in model", 
				omxQuotes(flatModel@name), "does not contain identical",
				"row and column dimnames.")
			stop(msg, call.=FALSE)
	}
	if (!isS4(means) && is.na(means)) return()
	means <- dimnames(means)
	if (is.null(means)) {
			msg <- paste("The expected means matrix associated",
				"with the FIML objective function in model", 
				omxQuotes(flatModel@name), "does not contain dimnames.")
			stop(msg, call.=FALSE)	
	}
	meanRows <- means[[1]]
	meanCols <- means[[2]]
	if (!is.null(meanRows) && length(meanRows) > 1) {
		msg <- paste("The expected means matrix associated",
			"with the FIML objective in model", 
			omxQuotes(flatModel@name), "is not a 1 x N matrix.")
			stop(msg, call.=FALSE)
	}
	if ((length(covCols) != length(meanCols)) || !all(covCols == meanCols)) {
			msg <- paste("The expected covariance and expected",
				"means matrices associated",
				"with the FIML objective function in model", 
				omxQuotes(flatModel@name), "do not contain identical",
				"dimnames.")
			stop(msg, call.=FALSE)
	}
}

generateDataRow <- function(flatModel, covName, dataName) {
	retval <- c()
	definitionNames <- dimnames(flatModel@datasets[[dataName]]@matrix)[[2]]
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


mxFIMLObjective <- function(covariance, means) {
	if (missing(covariance) || typeof(covariance) != "character") {
		stop("Covariance argument is not a string (the name of the expected covariance matrix)")
	}
	if (missing(means) || typeof(means) != "character") {
		stop("Means argument is not a string (the name of the expected means vector)")
	}
	return(new("MxFIMLObjective", covariance, means))
}

displayFIMLObjective <- function(objective) {
	cat("MxFIMLObjective", omxQuotes(objective@name), '\n')
	cat("covariance :", omxQuotes(objective@covariance), '\n')
	cat("means :", omxQuotes(objective@means), '\n')
	if (length(objective@result) == 0) {
		cat("Result : empty\n")
	} else {
		cat("Result : \n") 
		print(objective@result)
	}
	invisible(objective)
}


setMethod("print", "MxFIMLObjective", function(x,...) { 
	displayFIMLObjective(x) 
})

setMethod("show", "MxFIMLObjective", function(object) { 
	displayFIMLObjective(object) 
})
