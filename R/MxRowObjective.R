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


setClass(Class = "MxRowObjective",
	representation = representation(
		rowAlgebra = "MxCharOrNumber",
		rowResults = "MxCharOrNumber",
		filteredDataRow = "MxCharOrNumber",
		existenceVector = "MxCharOrNumber",
		reduceAlgebra = "MxCharOrNumber",
		definitionVars = "list",
		dims = "character",
		dataColumns = "numeric"),
	contains = "MxBaseObjective")

setMethod("initialize", "MxRowObjective",
	function(.Object, rowAlgebra, rowResults, filteredDataRow, 
		existenceVector, reduceAlgebra, dims,
		data = as.integer(NA), definitionVars = list(), name = 'objective') {
		.Object@name <- name
		.Object@rowAlgebra <- rowAlgebra
		.Object@rowResults <- rowResults
		.Object@reduceAlgebra <- reduceAlgebra
		.Object@filteredDataRow <- filteredDataRow
		.Object@existenceVector <- existenceVector
		.Object@data <- data
		.Object@definitionVars <- definitionVars
		.Object@dims <- dims
		return(.Object)
	}
)

setMethod("genericObjNewEntities", signature("MxRowObjective"),
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

setMethod("genericObjDependencies", signature("MxRowObjective"),
	function(.Object, dependencies) {
	sources <- c(.Object@reduceAlgebra)
	sources <- sources[!is.na(sources)]
	if (length(sources) > 0) {
		dependencies <- imxAddDependency(sources, .Object@name, dependencies)
	}
	return(dependencies)
})


setMethod("genericObjFunNamespace", signature("MxRowObjective"), 
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

setMethod("genericObjRename", signature("MxRowObjective"),
	function(.Object, oldname, newname) {
		.Object@rowAlgebra <- renameReference(.Object@rowAlgebra, oldname, newname)
		.Object@reduceAlgebra <- renameReference(.Object@reduceAlgebra, oldname, newname)
		.Object@data <- renameReference(.Object@data, oldname, newname)
		return(.Object)
})

setMethod("genericObjFunConvert", signature("MxRowObjective"), 
	function(.Object, flatModel, model, defVars) {
		modelname <- imxReverseIdentifier(model, .Object@name)[[1]]
		name <- .Object@name
		dataName <- .Object@data 
		if(is.na(dataName)) {
			msg <- paste("The MxRowObjective objective function",
				"does not have a dataset associated with it in model",
				omxQuotes(modelname))
			stop(msg, call. = FALSE)
		}
		mxDataObject <- flatModel@datasets[[dataName]]
		if (mxDataObject@type != 'raw') {
			msg <- paste("The dataset associated with the MxRowObjective objective", 
				"in model", omxQuotes(modelname), "is not raw data.")
			stop(msg, call. = FALSE)
		}
		.Object@definitionVars <- imxFilterDefinitionVariables(defVars, dataName)
		.Object@rowAlgebra <- imxLocateIndex(flatModel, .Object@rowAlgebra, name)
		.Object@rowResults <- imxLocateIndex(flatModel, .Object@rowResults, name)
		.Object@filteredDataRow <- imxLocateIndex(flatModel, .Object@filteredDataRow, name)
		.Object@existenceVector <- imxLocateIndex(flatModel, .Object@existenceVector, name)
		.Object@reduceAlgebra <- imxLocateIndex(flatModel, .Object@reduceAlgebra, name)
		.Object@data <- imxLocateIndex(flatModel, dataName, name)
		.Object@dataColumns <- generateRowDataColumns(flatModel, .Object@dims, dataName)
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
				"of the FIML objective function in model",
				omxQuotes(flatModel@name),
				"cannot be found in the column names of the data.")
			stop(msg, call. = FALSE)
		}
		retval[[i]] <- index - 1
	}
	return(retval)
}

setMethod("genericObjModelConvert", "MxRowObjective",
	function(.Object, job, model, namespace, flatJob) {
		rowAlgebraName <- .Object@rowAlgebra
		rowResultsName <- .Object@rowResults
		filteredDataRowName <- .Object@filteredDataRow
		existenceVectorName <- .Object@existenceVector
		reduceAlgebraName <- .Object@reduceAlgebra
		dimnames <- .Object@dims

		# Create the filtered data row
		if (is.na(filteredDataRowName)) {
			filteredDataRowName <- imxIdentifier(model@name, imxUntitledName())
		}
		filteredDataRow <- job[[filteredDataRowName]]
		if (!is.null(filteredDataRow)) {
			msg <- paste("The filteredDataRow cannot have name", 
				omxQuotes(simplifyName(filteredDataRowName, model@name)), 
				"because this entity already exists in the model")
			stop(msg, call. = FALSE)
		}
		filteredDataRow <- mxMatrix('Full', nrow = 1, ncol = length(dimnames))
		job[[filteredDataRowName]] <- filteredDataRow
		flatJob[[filteredDataRowName]] <- filteredDataRow
		pair <- imxReverseIdentifier(model, filteredDataRowName)
		if (model@name == pair[[1]]) {
			job[[.Object@name]]@filteredDataRow <- pair[[2]]
		} else {
			job[[.Object@name]]@filteredDataRow <- filteredDataRowName
		}

		# Create the existence vector
		if (is.na(existenceVectorName)) {
			existenceVectorName <- imxIdentifier(model@name, imxUntitledName())
		}
		existenceVector <- job[[existenceVectorName]]
		if (!is.null(existenceVector)) {
			msg <- paste("The existenceVector cannot have name", 
				omxQuotes(simplifyName(existenceVectorName, model@name)), 
				"because this entity already exists in the model")
			stop(msg, call. = FALSE)
		}
		existenceVector <- mxMatrix('Full', nrow = 1, ncol = length(dimnames), values = 1)
		job[[existenceVectorName]] <- existenceVector
		flatJob[[existenceVectorName]] <- existenceVector
		pair <- imxReverseIdentifier(model, existenceVectorName)
		if (model@name == pair[[1]]) {
			job[[.Object@name]]@existenceVector <- pair[[2]]
		} else {
			job[[.Object@name]]@existenceVector <- existenceVectorName
		}

		# Locate the row algebra
		rowAlgebra <- job[[rowAlgebraName]]
		if (is.null(rowAlgebra)) {
			msg <- paste("The rowAlgebra with name", 
				omxQuotes(simplifyName(rowAlgebraName, model@name)), 
				"is not defined in the model")
			stop(msg, call. = FALSE)
		}
		labelsData <- imxGenerateLabels(job)
		result <- evaluateMxObject(rowAlgebraName, flatJob, labelsData)
		if (nrow(result) != 1) {
			msg <- paste("The rowAlgebra with name", 
				omxQuotes(simplifyName(rowAlgebraName, model@name)), 
				"does not evaluate to a row vector")
			stop(msg, call. = FALSE)			
		}
		if (is.na(.Object@data)) {
			msg <- paste("The MxRowObjective objective function",
				"does not have a dataset associated with it in model",
				omxQuotes(model@name))
			stop(msg, call.=FALSE)		
		}
		mxDataObject <- flatJob@datasets[[.Object@data]]

		# Create the row results
		rows <- nrow(mxDataObject@observed)
		cols <- ncol(result)
		if (is.na(rowResultsName)) {
			rowResultsName <- imxIdentifier(model@name, imxUntitledName())
		}
		rowResults <- job[[rowResultsName]]
		if (!is.null(rowResults)) {
			msg <- paste("The rowResults cannot have name", 
				omxQuotes(simplifyName(rowResultsName, model@name)), 
				"because this entity already exists in the model")
			stop(msg, call. = FALSE)
		}
		rowResults <- mxMatrix('Full', nrow = rows, ncol = cols)
		job[[rowResultsName]] <- rowResults
		flatJob[[rowResultsName]] <- rowResults
		pair <- imxReverseIdentifier(model, rowResultsName)
		if (model@name == pair[[1]]) {
			job[[.Object@name]]@rowResults <- pair[[2]]
		} else {
			job[[.Object@name]]@rowResults <- rowResultsName
		}

		# Locate the reduce algebra
		if (is.na(reduceAlgebraName)) {
			reduceAlgebraName <- imxUntitledName()
			reduceAlgebra <- eval(substitute(mxAlgebra(x, reduceAlgebraName), 
				list(x = quote(as.symbol(rowResultsName)))))
			job[[model@name]][[reduceAlgebraName]] <- reduceAlgebra
			job[[.Object@name]]@reduceAlgebra <- reduceAlgebraName
			reduceAlgebraName <- imxIdentifier(model@name, reduceAlgebraName)
		}
		reduceAlgebra <- job[[reduceAlgebraName]]
		if (is.null(reduceAlgebra)) {
			msg <- paste("The reduceAlgebra with name", 
				omxQuotes(simplifyName(reduceAlgebraName, model@name)), 
				"is not defined in the model")
			stop(msg, call. = FALSE)
		}
		job@.newobjects <- TRUE
		job@.newobjective <- FALSE
		job@.newtree <- FALSE
		return(job)
	}
)

setMethod("genericObjInitialMatrix", "MxRowObjective",
	function(.Object, flatModel) {
		reduceAlgebraName <- .Object@reduceAlgebra
		labelsData <- imxGenerateLabels(flatModel)
		result <- evaluateMxObject(reduceAlgebraName, flatModel, labelsData)
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

mxRowObjective <- function(rowAlgebra, dimnames, rowResults = NA, filteredDataRow = NA, existenceVector = NA, reduceAlgebra = NA) {
	if (missing(rowAlgebra) || typeof(rowAlgebra) != "character") {
		stop("the 'rowAlgebra' argument is not a string (the name of the row-by-row algebra)")
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
	if (single.na(reduceAlgebra)) {
		reduceAlgebra <- as.character(NA)
	}
	if (!(is.vector(reduceAlgebra) && 
		(typeof(reduceAlgebra) == 'character') && 
		(length(reduceAlgebra) == 1))) {
		stop("the 'reduceAlgebra' argument is not a string (the name of the reduction algebra)")
	}
	return(new("MxRowObjective", rowAlgebra, rowResults, filteredDataRow, existenceVector, reduceAlgebra, dimnames))
}

printSlot <- function(object, slotName) {
	val <- slot(object, slotName)
	if (single.na(val)) {
		msg <- paste('@', slotName, ' : NA \n', sep = '')
	} else {
		msg <- paste('@', slotName, ' : ',omxQuotes(val), '\n', sep = '')
	}
	cat(msg)
}

displayRowObjective <- function(objective) {
	cat("MxRowObjective", omxQuotes(objective@name), '\n')
	cat("@rowAlgebra :", omxQuotes(objective@rowAlgebra), '\n')
	printSlot(objective, "rowResults")
	printSlot(objective, "filteredDataRow")
	printSlot(objective, "existenceVector")
	printSlot(objective, "reduceAlgebra")
	if (length(objective@result) == 0) {
		cat("@result: (not yet computed) ")
	} else {
		cat("@result:\n")
	}
	print(objective@result)
	invisible(objective)
}


setMethod("print", "MxRowObjective", function(x,...) { 
	displayRowObjective(x) 
})

setMethod("show", "MxRowObjective", function(object) { 
	displayRowObjective(object) 
})
