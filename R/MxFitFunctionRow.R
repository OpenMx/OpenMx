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


setClass(Class = "MxFitFunctionRow",
	representation = representation(
		rowAlgebra = "MxCharOrNumber",
		rowResults = "MxCharOrNumber",
	    units = "character",
		filteredDataRow = "MxCharOrNumber",
		existenceVector = "MxCharOrNumber",
		reduceAlgebra = "MxCharOrNumber",
		data = "MxCharOrNumber",
		dims = "character",
		dataColumns = "numeric",
		dataRowDeps = "integer"),
	contains = "MxBaseFitFunction")

setMethod("initialize", "MxFitFunctionRow",
	function(.Object, rowAlgebra, rowResults, units, filteredDataRow, 
		existenceVector, reduceAlgebra, dims,
		data = as.integer(NA), name = 'fitfunction') {
		.Object@name <- name
		.Object@rowAlgebra <- rowAlgebra
		.Object@rowResults <- rowResults
		.Object@units <- units
		.Object@reduceAlgebra <- reduceAlgebra
		.Object@filteredDataRow <- filteredDataRow
		.Object@existenceVector <- existenceVector
		.Object@data <- data
		.Object@dims <- dims
		.Object@expectation <- as.integer(NA)
		return(.Object)
	}
)

setMethod("genericFitNewEntities", signature("MxFitFunctionRow"),
	function(.Object) {
		if (is.na(.Object@rowResults) && is.na(.Object@filteredDataRow) && is.na(.Object@existenceVector)) {
			return(NULL)
		} else {
			a <- .Object@rowResults
			b <- .Object@filteredDataRow
			c <- .Object@existenceVector
			retval <- c(a, b, c)
			retval <- as.character(na.omit(retval))
			return(retval)
		}
	}
)

setMethod("genericFitDependencies", signature("MxFitFunctionRow"),
	function(.Object, flatModel, dependencies) {
		dependencies <- callNextMethod()
		reduceAlgebra <- .Object@reduceAlgebra
		rowAlgebra <- .Object@rowAlgebra
		rowResults <- .Object@rowResults
		dependencies <- imxAddDependency(reduceAlgebra, .Object@name, dependencies)
		dependencies <- imxAddDependency(rowAlgebra, rowResults, dependencies)		
		return(dependencies)
})


setMethod("qualifyNames", signature("MxFitFunctionRow"), 
	function(.Object, modelname, namespace) {
		.Object@name <- imxIdentifier(modelname, .Object@name)
		.Object@rowAlgebra <- imxConvertIdentifier(.Object@rowAlgebra, 
			modelname, namespace)
		.Object@rowResults <- imxConvertIdentifier(.Object@rowResults,
			modelname, namespace)
		.Object@filteredDataRow <- imxConvertIdentifier(.Object@filteredDataRow, 
			modelname, namespace)
		.Object@existenceVector <- imxConvertIdentifier(.Object@existenceVector, 
			modelname, namespace)
		.Object@reduceAlgebra <- imxConvertIdentifier(.Object@reduceAlgebra, 
			modelname, namespace)
		.Object@data <- imxConvertIdentifier(.Object@data, 
			modelname, namespace)
		return(.Object)
})

setMethod("genericFitRename", signature("MxFitFunctionRow"),
	function(.Object, oldname, newname) {
		.Object@rowAlgebra <- renameReference(.Object@rowAlgebra, oldname, newname)
		.Object@reduceAlgebra <- renameReference(.Object@reduceAlgebra, oldname, newname)
		.Object@data <- renameReference(.Object@data, oldname, newname)
		return(.Object)
})

setMethod("genericFitFunConvert", signature("MxFitFunctionRow"), 
	function(.Object, flatModel, model, labelsData, dependencies) {
		modelname <- imxReverseIdentifier(model, .Object@name)[[1]]
		name <- .Object@name
		dataName <- .Object@data
		if(is.na(dataName)) {
			msg <- paste("The MxFitFunctionRow fit function",
				"does not have a dataset associated with it in model",
				omxQuotes(modelname))
			stop(msg, call. = FALSE)
		}
		mxDataObject <- flatModel@datasets[[dataName]]
		if (mxDataObject@type != 'raw') {
			msg <- paste("The dataset associated with the MxFitFunctionRow fit function", 
				"in model", omxQuotes(modelname), "is not raw data.")
			stop(msg, call. = FALSE)
		}
		dataRowDeps <- union(dependencies[[.Object@filteredDataRow]], dependencies[[.Object@existenceVector]])
		dataRowDeps <- sapply(dataRowDeps, doLocateIndex, flatModel, flatModel@name, USE.NAMES=FALSE)
		dataRowDeps <- as.integer(dataRowDeps)
		.Object@rowAlgebra <- imxLocateIndex(flatModel, .Object@rowAlgebra, name)
		.Object@rowResults <- imxLocateIndex(flatModel, .Object@rowResults, name)
		.Object@filteredDataRow <- imxLocateIndex(flatModel, .Object@filteredDataRow, name)
		.Object@existenceVector <- imxLocateIndex(flatModel, .Object@existenceVector, name)
		.Object@reduceAlgebra <- imxLocateIndex(flatModel, .Object@reduceAlgebra, name)
		.Object@data <- imxLocateIndex(flatModel, dataName, name)
		.Object@dataColumns <- generateRowDataColumns(flatModel, .Object@dims, dataName)
		.Object@dataRowDeps <- dataRowDeps
		if (length(mxDataObject@observed) == 0) {
			.Object@data <- as.integer(NA)
		}
		return(.Object)
})

generateRowDataColumns <- function(flatModel, expectedNames, dataName) {
	retval <- c()
	dataColumnNames <- dimnames(flatModel@datasets[[dataName]]@observed)[[2]]
	for(i in 1:length(expectedNames)) {
		targetName <- expectedNames[[i]]
		index <- match(targetName, dataColumnNames)
		if(is.na(index)) {
			msg <- paste("The column name", omxQuotes(targetName),
				"in the expected covariance matrix",
				"of the MxFitFunctionRow fit function in model",
				omxQuotes(flatModel@name),
				"cannot be found in the column names of the data.")
			stop(msg, call. = FALSE)
		}
		retval[[i]] <- index - 1
	}
	return(retval)
}

setMethod("genericFitAddEntities", "MxFitFunctionRow",
	function(.Object, job, flatJob, labelsData) {
		name <- .Object@name
		modelname <- job@name
		rowAlgebraName <- .Object@rowAlgebra
		rowResultsName <- .Object@rowResults
		filteredDataRowName <- .Object@filteredDataRow
		existenceVectorName <- .Object@existenceVector
		reduceAlgebraName <- .Object@reduceAlgebra
		dimnames <- .Object@dims
		job@.newobjects <- TRUE

		# Create the filtered data row
		filteredDataRow <- flatJob[[filteredDataRowName]]
		if (!is.null(filteredDataRow)) {
			msg <- paste("The filteredDataRow cannot have name", 
				omxQuotes(filteredDataRowName), 
				"because this entity already exists in the model")
			stop(msg, call. = FALSE)
		}
		filteredDataRow <- mxMatrix('Full', nrow = 1, ncol = length(dimnames))
    filteredDataRow@.persist <- FALSE
		job[[filteredDataRowName]] <- filteredDataRow
		flatJob[[filteredDataRowName]] <- filteredDataRow

		# Create the existence vector
		if (!is.na(existenceVectorName)) {
			existenceVector <- job[[existenceVectorName]]
			if (!is.null(existenceVector)) {
				msg <- paste("The existenceVector cannot have name", 
					omxQuotes(existenceVectorName), 
					"because this entity already exists in the model")
				stop(msg, call. = FALSE)
			}
			existenceVector <- mxMatrix('Full', nrow = 1, ncol = length(dimnames), values = 1)
			existenceVector@.persist <- FALSE
      job[[existenceVectorName]] <- existenceVector
			flatJob[[existenceVectorName]] <- existenceVector
		}

		# Locate the row algebra
		rowAlgebra <- job[[rowAlgebraName]]
		if (is.null(rowAlgebra)) {
			msg <- paste("The rowAlgebra with name", 
				omxQuotes(rowAlgebraName), 
				"is not defined in the model")
			stop(msg, call. = FALSE)
		}
		tuple <- evaluateMxObject(rowAlgebraName, flatJob, labelsData, new.env(parent = emptyenv()))
		result <- tuple[[1]]
		if (nrow(result) != 1) {
			msg <- paste("The rowAlgebra with name", 
				omxQuotes(rowAlgebraName), 
				"does not evaluate to a row vector")
			stop(msg, call. = FALSE)			
		}
		if (is.na(.Object@data)) {
			msg <- paste("The MxFitFunctionRow fit function",
				"does not have a dataset associated with it in model",
				omxQuotes(modelname))
			stop(msg, call.=FALSE)		
		}
		mxDataObject <- flatJob@datasets[[.Object@data]]

		# Create the row results
		rows <- nrow(mxDataObject@observed)
		cols <- ncol(result)
		rowResults <- job[[rowResultsName]]
		if (!is.null(rowResults)) {
			msg <- paste("The rowResults cannot have name", 
				omxQuotes(rowResultsName), 
				"because this entity already exists in the model")
			stop(msg, call. = FALSE)
		}
		rowResults <- mxMatrix('Full', nrow = rows, ncol = cols)
		rowResults@.persist <- FALSE
    job[[rowResultsName]] <- rowResults

		# Locate the reduce algebra
		reduceAlgebra <- job[[reduceAlgebraName]]
		if (is.null(reduceAlgebra)) {
			msg <- paste("The reduceAlgebra with name", 
				omxQuotes(reduceAlgebraName), 
				"is not defined in the model")
			stop(msg, call. = FALSE)
		}
		return(job)
	}
)

setMethod("genericFitInitialMatrix", "MxFitFunctionRow",
	function(.Object, flatModel) {
		reduceAlgebraName <- .Object@reduceAlgebra
		labelsData <- imxGenerateLabels(flatModel)
		tuple <- evaluateMxObject(reduceAlgebraName, flatModel, labelsData, new.env(parent = emptyenv()))
		result <- tuple[[1]]
		return(result)
	}
)

checkStringArgument <- function(arg, name) {
	if (single.na(arg)) {
		arg <- as.character(NA)	
	} else if (length(unlist(strsplit(arg, imxSeparatorChar, fixed = TRUE))) > 1) {
		stop(paste("the", omxQuotes(name), "argument cannot contain the", 
			omxQuotes(imxSeparatorChar), 
			"character"))
	}
	if (!(is.vector(arg) && 
		(typeof(arg) == 'character') && 
		(length(arg) == 1))) {
		stop("the", omxQuotes(name), "argument is not a string")
	}
	return(arg)
}

mxFitFunctionRow <- function(rowAlgebra, reduceAlgebra, dimnames, rowResults = "rowResults", 
	filteredDataRow = "filteredDataRow", existenceVector = "existenceVector", units="-2lnL") {
	if (missing(rowAlgebra) || typeof(rowAlgebra) != "character") {
		stop("the 'rowAlgebra' argument is not a string (the name of the row-by-row algebra)")
	}
	if (missing(reduceAlgebra) || typeof(reduceAlgebra) != "character") {
		stop("the 'reduceAlgebra' argument is not a string (the name of the reduction algebra)")
	}
	if (missing(dimnames) || typeof(dimnames) != "character") {
		stop("the 'dimnames' argument is not a string (the column names from the data set)")
	}
	if (any(is.na(dimnames))) {
		stop("NA values are not allowed for 'dimnames' vector")
	}
	rowResults <- checkStringArgument(rowResults, "rowResults")
	filteredDataRow <- checkStringArgument(filteredDataRow, "filteredDataRow")
	existenceVector <- checkStringArgument(existenceVector, "existenceVector")
	return(new("MxFitFunctionRow", rowAlgebra, rowResults, units,
		   filteredDataRow, existenceVector, reduceAlgebra, dimnames))
}

printSlot <- function(object, slotName) {
	val <- slot(object, slotName)
	if (single.na(val)) {
		msg <- paste('$', slotName, ' : NA \n', sep = '')
	} else {
		msg <- paste('$', slotName, ' : ',omxQuotes(val), '\n', sep = '')
	}
	cat(msg)
}

displayRowFitFunction <- function(fitfunction) {
	cat("MxFitFunctionRow", omxQuotes(fitfunction@name), '\n')
	cat("$rowAlgebra :", omxQuotes(fitfunction@rowAlgebra), '\n')
	cat("$units: ", omxQuotes(fitfunction@units), '\n')
	printSlot(fitfunction, "rowResults")
	printSlot(fitfunction, "filteredDataRow")
	printSlot(fitfunction, "existenceVector")
	printSlot(fitfunction, "reduceAlgebra")
	if (length(fitfunction@result) == 0) {
		cat("$result: (not yet computed) ")
	} else {
		cat("$result:\n")
	}
	print(fitfunction@result)
	invisible(fitfunction)
}


setMethod("print", "MxFitFunctionRow", function(x,...) { 
	displayRowFitFunction(x) 
})

setMethod("show", "MxFitFunctionRow", function(object) { 
	displayRowFitFunction(object) 
})
