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
		reduceAlgebra = "MxCharOrNumber",
		definitionVars = "list"),
	contains = "MxBaseObjective")

setMethod("initialize", "MxRowObjective",
	function(.Object, rowAlgebra, rowResults, reduceAlgebra,
		data = as.integer(NA), definitionVars = list(), name = 'objective') {
		.Object@name <- name
		.Object@rowAlgebra <- rowAlgebra
		.Object@rowResults <- rowResults
		.Object@reduceAlgebra <- reduceAlgebra		
		.Object@data <- data
		.Object@definitionVars <- definitionVars
		return(.Object)
	}
)

setMethod("omxObjDependencies", signature("MxRowObjective"),
	function(.Object, dependencies) {
	sources <- c(.Object@reduceAlgebra)
	sources <- sources[!is.na(sources)]
	if (length(sources) > 0) {
		dependencies <- omxAddDependency(sources, .Object@name, dependencies)
	}
	return(dependencies)
})


setMethod("omxObjFunNamespace", signature("MxRowObjective"), 
	function(.Object, modelname, namespace) {
		.Object@name <- omxIdentifier(modelname, .Object@name)
		.Object@rowAlgebra <- omxConvertIdentifier(.Object@rowAlgebra, 
			modelname, namespace)
		.Object@rowResults <- omxConvertIdentifier(.Object@rowResults, 
			modelname, namespace)
		.Object@reduceAlgebra <- omxConvertIdentifier(.Object@reduceAlgebra, 
			modelname, namespace)
		.Object@data <- omxConvertIdentifier(.Object@data, 
			modelname, namespace)
		return(.Object)
})

setMethod("omxObjFunConvert", signature("MxRowObjective"), 
	function(.Object, flatModel, model) {
		modelname <- omxReverseIdentifier(model, .Object@name)[[1]]
		name <- .Object@name
		if(is.na(.Object@data)) {
			msg <- paste("The MxRowObjective objective function",
				"does not have a dataset associated with it in model",
				omxQuotes(modelname))
			stop(msg, call.=FALSE)
		}
		mxDataObject <- flatModel@datasets[[.Object@data]]
		if (mxDataObject@type != 'raw') {
			msg <- paste("The dataset associated with the MxRowObjective objective", 
				"in model", omxQuotes(modelname), "is not raw data.")
			stop(msg, call.=FALSE)
		}
		.Object@rowAlgebra <- omxLocateIndex(flatModel, .Object@rowAlgebra, name)
		.Object@rowResults <- omxLocateIndex(flatModel, .Object@rowResults, name)
		.Object@reduceAlgebra <- omxLocateIndex(flatModel, .Object@reduceAlgebra, name)
		.Object@definitionVars <- generateDefinitionList(flatModel)
		if (length(mxDataObject@observed) == 0) {
			.Object@data <- as.integer(NA)
		}
		return(.Object)
})

setMethod("omxObjModelConvert", "MxRowObjective",
	function(.Object, job, model, flatJob) {
		rowAlgebraName <- .Object@rowAlgebra
		rowResultsName <- .Object@rowResults
		reduceAlgebraName <- .Object@reduceAlgebra
		rowAlgebra <- job[[rowAlgebraName]]
		if (is.null(rowAlgebra)) {
			msg <- paste("The rowAlgebra with name", 
				omxQuotes(simplifyName(rowAlgebraName, model@name)), 
				"is not defined in the model")
			stop(msg, call. = FALSE)
		}
		labelsData <- omxGenerateLabels(job)		
		result <- tryCatch(eval(computeSymbol(as.symbol(rowAlgebraName), flatJob, labelsData)), 
			error = function(x) {
				stop(paste("The entity", 
					omxQuotes(simplifyName(rowAlgebraName, model@name)), 
					"in model", omxQuotes(model@name), 
					"generated the error message:",
					x$message), call. = FALSE)
		})
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
		rows <- nrow(mxDataObject@observed)
		cols <- ncol(result)
		if (is.na(rowResultsName)) {
			rowResultsName <- omxUntitledName()			
		}
		rowResults <- job[[model@name]][[rowResultsName]]
		if (!is.null(rowResults)) {
			msg <- paste("The rowResults cannot have name", 
				omxQuotes(simplifyName(rowAlgebraName, model@name)), 
				"because this entity already exists in the model")
			stop(msg, call. = FALSE)
		}
		rowResults <- mxMatrix('Full', nrow = rows, ncol = cols, name = rowResultsName)
		job[[model@name]][[rowResultsName]] <- rowResults
		job[[.Object@name]]@rowResults <- rowResultsName
		if (is.na(reduceAlgebraName)) {
			reduceAlgebraName <- omxUntitledName()
			reduceAlgebra <- eval(substitute(mxAlgebra(x, reduceAlgebraName), 
				list(x = as.symbol(rowResultsName))))
			job[[model@name]][[reduceAlgebraName]] <- reduceAlgebra		
			job[[.Object@name]]@reduceAlgebra <- reduceAlgebraName
			reduceAlgebraName <- omxIdentifier(model@name, reduceAlgebraName)
		}
		reduceAlgebra <- job[[reduceAlgebraName]]
		if (is.null(reduceAlgebra)) {
			msg <- paste("The reduceAlgebra with name", 
				omxQuotes(simplifyName(reduceAlgebraName, model@name)), 
				"is not defined in the model")
			stop(msg, call. = FALSE)
		}
		return(job)
	}
)

setMethod("omxObjInitialMatrix", "MxRowObjective",
	function(.Object, flatModel) {
		reduceAlgebraName <- .Object@reduceAlgebra
		labelsData <- omxGenerateLabels(flatModel)		
		result <- tryCatch(eval(computeSymbol(as.symbol(reduceAlgebraName), flatModel, labelsData)), 
			error = function(x) {
				stop(paste("The entity", 
					omxQuotes(simplifyName(reduceAlgebraName, flatModel@name)), 
					"in model", omxQuotes(flatModel@name), 
					"generated the error message:",
					x$message), call. = FALSE)
		})
		return(result)
	}
)

mxRowObjective <- function(rowAlgebra, rowResults = NA, reduceAlgebra = NA) {
	if (missing(rowAlgebra) || typeof(rowAlgebra) != "character") {
		stop("the 'rowAlgebra' argument is not a string (the name of the row-by-row algebra)")
	}
	if (single.na(rowResults)) {
		rowResults <- as.character(NA)	
	} else if (length(unlist(strsplit(rowResults, omxSeparatorChar, fixed = TRUE))) > 1) {
		stop(paste("the 'rowResults' argument cannot contain the", 
			omxQuotes(omxSeparatorChar), 
			"character"))
	}
	if (!(is.vector(rowResults) && 
		(typeof(rowResults) == 'character') && 
		(length(rowResults) == 1))) {
		stop("the 'rowResults' argument is not a string (the name for the results matrix)")
	}
	if (single.na(reduceAlgebra)) {
		reduceAlgebra <- as.character(NA)	
	}
	if (!(is.vector(reduceAlgebra) && 
		(typeof(reduceAlgebra) == 'character') && 
		(length(reduceAlgebra) == 1))) {
		stop("the 'reduceAlgebra' argument is not a string (the name of the reduction algebra)")
	}
	return(new("MxRowObjective", rowAlgebra, rowResults, reduceAlgebra))
}

displayRowObjective <- function(objective) {
	cat("MxRowObjective", omxQuotes(objective@name), '\n')
	cat("@rowAlgebra :", omxQuotes(objective@rowAlgebra), '\n')
	if (single.na(objective@rowResults)) {
		cat("@rowResults : NA \n")
	} else {
		cat("@rowResults :", omxQuotes(objective@rowResults), '\n')	
	}
	if (single.na(objective@reduceAlgebra)) {
		cat("@reduceAlgebra : NA \n")
	} else {
		cat("@reduceAlgebra :", omxQuotes(objective@reduceAlgebra), '\n')	
	}	
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
