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


setClass(Class = "MxExpectationRAM",
	representation = representation(
		A = "MxCharOrNumber",
		S = "MxCharOrNumber",
		F = "MxCharOrNumber",
		M = "MxCharOrNumber",
		thresholds = "MxCharOrNumber",
		dims = "character",
		definitionVars = "list",
		dataColumns = "numeric",
		thresholdColumns = "numeric",
		thresholdLevels = "numeric",
		depth = "integer",
		threshnames = "character",
		usePPML = "logical",
		ppmlData = "MxData",
		UnfilteredExpCov = "matrix",
	        numStats = "numeric"),
	contains = "MxBaseExpectation")

setMethod("initialize", "MxExpectationRAM",
	function(.Object, A, S, F, M, dims, thresholds, threshnames,
		data = as.integer(NA), name = 'expectation') {
		.Object@name <- name
		.Object@A <- A
		.Object@S <- S
		.Object@F <- F
		.Object@M <- M
		.Object@data <- data
		.Object@dims <- dims
		.Object@thresholds <- thresholds
		.Object@definitionVars <- list()
		.Object@threshnames <- threshnames
		.Object@usePPML <- FALSE
		.Object@UnfilteredExpCov <- matrix()
		return(.Object)
	}
)

setMethod("genericExpDependencies", signature("MxExpectationRAM"),
	function(.Object, dependencies) {
	sources <- c(.Object@A, .Object@S, .Object@F, .Object@M, .Object@thresholds)
	sources <- sources[!is.na(sources)]
	dependencies <- imxAddDependency(sources, .Object@name, dependencies)
	return(dependencies)
})

setMethod("qualifyNames", signature("MxExpectationRAM"), 
	function(.Object, modelname, namespace) {
		.Object@name <- imxIdentifier(modelname, .Object@name)
		.Object@A <- imxConvertIdentifier(.Object@A, modelname, namespace)
		.Object@S <- imxConvertIdentifier(.Object@S, modelname, namespace)
		.Object@F <- imxConvertIdentifier(.Object@F, modelname, namespace)
		.Object@M <- imxConvertIdentifier(.Object@M, modelname, namespace)
		.Object@data <- imxConvertIdentifier(.Object@data, modelname, namespace)
		.Object@thresholds <- sapply(.Object@thresholds, 
			imxConvertIdentifier, modelname, namespace)
		return(.Object)
})

setMethod("genericExpRename", signature("MxExpectationRAM"),
	function(.Object, oldname, newname) {
		.Object@A <- renameReference(.Object@A, oldname, newname)
		.Object@S <- renameReference(.Object@S, oldname, newname)
		.Object@F <- renameReference(.Object@F, oldname, newname)
		.Object@M <- renameReference(.Object@M, oldname, newname)
		.Object@data <- renameReference(.Object@data, oldname, newname)
		.Object@thresholds <- sapply(.Object@thresholds, renameReference, oldname, newname)		
		return(.Object)
})

setMethod("genericExpFunConvert", signature("MxExpectationRAM"), 
	function(.Object, flatModel, model, labelsData, defVars, dependencies) {
		modelname <- imxReverseIdentifier(model, .Object@name)[[1]]	
		name <- .Object@name
		aMatrix <- .Object@A
		sMatrix <- .Object@S
		fMatrix <- .Object@F
		mMatrix <- .Object@M
		data <- .Object@data
		if(is.na(data)) {
			msg <- paste("The RAM expectation function",
				"does not have a dataset associated with it in model",
				omxQuotes(modelname),
				"\nSee ?mxData() to see how to add data to your model")
			stop(msg, call. = FALSE)
		}
		mxDataObject <- flatModel@datasets[[.Object@data]]
		if(!is.na(mMatrix) && single.na(mxDataObject@means) && mxDataObject@type != "raw") {
			msg <- paste("The RAM expectation function",
				"has an expected means vector but",
				"no observed means vector in model",
				omxQuotes(modelname))
			stop(msg, call. = FALSE)
		}
		if(!single.na(mxDataObject@means) && is.null(flatModel[[mMatrix]])) {
			msg <- paste("The RAM expectation function",
				"has an observed means vector but",
				"no expected means vector in model",
				omxQuotes(modelname))
			stop(msg, call. = FALSE)		
		}
		checkNumericData(mxDataObject)
		.Object@A <- imxLocateIndex(flatModel, aMatrix, name)
		.Object@S <- imxLocateIndex(flatModel, sMatrix, name)
		.Object@F <- imxLocateIndex(flatModel, fMatrix, name)
		.Object@M <- imxLocateIndex(flatModel, mMatrix, name)
		.Object@data <- as.integer(imxLocateIndex(flatModel, data, name))
		verifyObservedNames(mxDataObject@observed, mxDataObject@means, mxDataObject@type, flatModel, modelname, "RAM")
		fMatrix <- flatModel[[fMatrix]]@values
		if (is.null(dimnames(fMatrix))) {
			msg <- paste("The F matrix of model",
				omxQuotes(modelname), "does not contain dimnames")
			stop(msg, call. = FALSE)
		}
		if (is.null(dimnames(fMatrix)[[2]])) {
			msg <- paste("The F matrix of model",
				omxQuotes(modelname), "does not contain colnames")
			stop(msg, call. = FALSE)
		}
		mMatrix <- flatModel[[mMatrix]]		
		if(!is.null(mMatrix)) {
			means <- dimnames(mMatrix)
			if (is.null(means)) {
				msg <- paste("The M matrix associated",
				"with the RAM expectation function in model", 
				omxQuotes(modelname), "does not contain dimnames.")
				stop(msg, call. = FALSE)	
			}
			meanRows <- means[[1]]
			meanCols <- means[[2]]
			if (!is.null(meanRows) && length(meanRows) > 1) {
				msg <- paste("The M matrix associated",
				"with the RAM expectation function in model", 
				omxQuotes(modelname), "is not a 1 x N matrix.")
				stop(msg, call. = FALSE)
			}
			if (!identical(dimnames(fMatrix)[[2]], meanCols)) {
				msg <- paste("The column names of the F matrix",
					"and the column names of the M matrix",
					"in model", 
					omxQuotes(modelname), "do not contain identical",
					"names.")
				stop(msg, call. = FALSE)
			}
		}
		translatedNames <- fMatrixTranslateNames(fMatrix, modelname)
		.Object@depth <- generateRAMDepth(flatModel, aMatrix, model@options)
		if (mxDataObject@type == 'raw' || mxDataObject@type == 'acov') {
			threshName <- .Object@thresholds
			checkNumberOrdinalColumns(mxDataObject)
			.Object@definitionVars <- imxFilterDefinitionVariables(defVars, data)
			.Object@dataColumns <- generateDataColumns(flatModel, translatedNames, data)
			verifyThresholds(flatModel, model, labelsData, data, translatedNames, threshName)
			.Object@thresholds <- imxLocateIndex(flatModel, threshName, name)
			retval <- generateThresholdColumns(flatModel, model, labelsData, translatedNames, data, threshName)
			.Object@thresholdColumns <- retval[[1]]
			.Object@thresholdLevels <- retval[[2]]
			if (length(mxDataObject@observed) == 0) {
				.Object@data <- as.integer(NA)
			}
			if (single.na(.Object@dims)) {
				.Object@dims <- translatedNames
			}
		} else {
			.Object@thresholds <- as.integer(NA)
			targetNames <- rownames(mxDataObject@observed)
			if (!identical(translatedNames, targetNames)) {
				varsNotInData <- translatedNames[!(translatedNames %in% targetNames)]
				msg <- paste("The names of the manifest",
					"variables in the F matrix of model",
					omxQuotes(modelname), "does not match the",
					"dimnames of the observed covariance matrix.")
				if (length(varsNotInData) > 0) {
					msg <- paste(msg,
						"To get you started, the following variables are used but",
						"are not in the observed data:", 
						omxQuotes(varsNotInData))
				}
				stop(msg, call. = FALSE)
			}
		}
		return(.Object)
})

generateRAMDepth <- function(flatModel, aMatrixName, modeloptions) {
	mxObject <- flatModel[[aMatrixName]]
	if (!is(mxObject, "MxMatrix")) {
		return(as.integer(NA))
	}
	if (identical(modeloptions[['RAM Inverse Optimization']], "No")) {
		return(as.integer(NA))
	}
	if (is.null(modeloptions[['RAM Inverse Optimization']]) &&
		identical(getOption('mxOptions')[['RAM Inverse Optimization']], "No")) {
		return(as.integer(NA))
	}	
	maxdepth <- modeloptions[['RAM Max Depth']]
	if (is.null(maxdepth) || (length(maxdepth) != 1) ||
		is.na(maxdepth) || !is.numeric(maxdepth) || maxdepth < 0) {
		maxdepth <- nrow(mxObject) - 1
	}
	return(omxGetRAMDepth(mxObject, maxdepth))
}

omxGetRAMDepth <- function(A, maxdepth = nrow(A) - 1) {
	mxObject <- A
	aValues <- matrix(0, nrow(mxObject), ncol(mxObject))
	defvars <- apply(mxObject@labels, c(1,2), imxIsDefinitionVariable)
	squarebrackets <- mxObject@.squareBrackets
	aValues[mxObject@free] <- 1
	aValues[mxObject@values != 0] <- 1
	aValues[defvars] <- 1
	aValues[squarebrackets] <- 1
	return(generateDepthHelper(aValues, aValues, 0, maxdepth))
}

generateDepthHelper <- function(aValues, currentProduct, depth, maxdepth) {
	if (depth > maxdepth) {
		return(as.integer(NA))
	}
	if (all(currentProduct == 0)) { 
		return(as.integer(depth))
	} else {
		return(generateDepthHelper(aValues, currentProduct %*% aValues, depth + 1, maxdepth))
	}
}

fMatrixTranslateNames <- function(fMatrix, modelName) {
	retval <- character()
	colNames <- dimnames(fMatrix)[[2]]
	for(i in 1:nrow(fMatrix)) {
		irow <- fMatrix[i,]
		matches <- which(irow == 1)
		if (length(matches) != 1) {
			err <- paste("The model",
				omxQuotes(modelName), "does not contain",
				"a valid F matrix")
			stop(err, call. = FALSE)
		}
		retval[[i]] <- colNames[[matches[[1]]]]
	}
	return(retval)
}

updateRAMdimnames <- function(flatExpectation, flatJob) {
	fMatrixName <- flatExpectation@F
	mMatrixName <- flatExpectation@M
	if (is.na(mMatrixName)) {
		mMatrix <- NA
	} else {
		mMatrix <- flatJob[[mMatrixName]]
	}
	fMatrix <- flatJob[[fMatrixName]]
	if (is.null(fMatrix)) {
		modelname <- getModelName(flatExpectation)
		stop(paste("Unknown F matrix name", 
			omxQuotes(simplifyName(fMatrixName, modelname)),
			"detected in the RAM expectation function",
			"of model", omxQuotes(modelname)), call. = FALSE)
	}
	dims <- flatExpectation@dims
	if (!is.null(dimnames(fMatrix)) && !single.na(dims) && 
		!identical(dimnames(fMatrix)[[2]], dims)) {
		modelname <- getModelName(flatExpectation)
		msg <- paste("The F matrix associated",
			"with the RAM expectation function in model", 
			omxQuotes(modelname), "contains dimnames and",
			"the expectation function has specified dimnames")
		stop(msg, call.=FALSE)		
	}
	if (is.null(dimnames(fMatrix)) && !single.na(dims)) {
		dimnames(flatJob[[fMatrixName]]) <- list(c(), dims)
	}

	if (!isS4(mMatrix) && (is.null(mMatrix) || is.na(mMatrix))) {
		return(flatJob)
	}

	if (!is.null(dimnames(mMatrix)) && !single.na(dims) &&
		!identical(dimnames(mMatrix), list(NULL, dims))) {
		modelname <- getModelName(flatExpectation)
		msg <- paste("The M matrix associated",
			"with the RAM expectation function in model", 
			omxQuotes(modelname), "contains dimnames and",
			"the expectation function has specified dimnames")
		stop(msg, call.=FALSE)	
	}

	if (is.null(dimnames(mMatrix)) && !single.na(dims)) {
		dimnames(flatJob[[mMatrixName]]) <- list(NULL, dims)
	}

	return(flatJob)
}

setMethod("genericExpAddEntities", "MxExpectationRAM",
	  function(.Object, job, flatJob, labelsData) {
		  precision <- "Function precision"
		  if(!single.na(.Object@thresholds)) {
			  if (is.null(job@options[[precision]])) {
				  job <- mxOption(job, precision, 1e-9)
			  }
		  }

		  ppmlModelOption <- job@options$UsePPML
		  if (is.null(ppmlModelOption)) {
			  enablePPML <- (getOption("mxOptions")$UsePPML == "Yes")
		  } else {
			  enablePPML <- (ppmlModelOption == "Yes")
		  }

		  if (enablePPML) {
			  aMatrix <- job[[.Object@A]]
			  aMatrixFixed <- !is.null(aMatrix) && is(aMatrix, "MxMatrix") && all(!aMatrix@free)
			  enablePPML <- aMatrixFixed
		  }

		  if (enablePPML) {
			  job <- PPMLTransformModel(job)
			  job@.newobjects <- TRUE
		  }

		  return(job)
	  })

setMethod("genericExpConvertEntities", "MxExpectationRAM",
	function(.Object, flatModel, namespace, labelsData) {
		if(is.na(.Object@data)) {
			modelname <- getModelName(.Object)
			msg <- paste("The RAM expectation function",
				"does not have a dataset associated with it in model",
				omxQuotes(modelname))
			stop(msg, call.=FALSE)
		}
		
		flatModel <- updateRAMdimnames(.Object, flatModel)

		if (flatModel@datasets[[.Object@data]]@type != 'raw') {
			return(flatModel)
		}

		if (is.na(.Object@M) || is.null(flatModel[[.Object@M]])) {
			modelname <- getModelName(.Object)
			msg <- paste("The RAM expectation function",
				"has raw data but is missing",
				"an expected means vector in model",
				omxQuotes(modelname))
			stop(msg, call.=FALSE)
		}
				
		flatModel <- updateThresholdDimnames(.Object, flatModel, labelsData)

		return(flatModel)
	}
)

##' imxSimpleRAMPredicate
##'
##' This is an internal function exported for those people who know
##' what they are doing.
##'
##' @param model model
imxSimpleRAMPredicate <- function(model) {
	if (is.null(model$expectation) || !is(model$expectation, "MxExpectationRAM")) {
		return(FALSE)
	}
	nameA <- model$expectation@A
	nameS <- model$expectation@S
	A <- model[[nameA]]
	S <- model[[nameS]]
	if (is.null(A) || is.null(S)) {
		return(FALSE)
	}
	return(is(A, "MxMatrix") && is(S, "MxMatrix"))
}

mxExpectationRAM <- function(A="A", S="S", F="F", M = NA, dimnames = NA, thresholds = NA,
	threshnames = dimnames) {

	if (typeof(A) != "character") {
		msg <- paste("argument 'A' is not a string",
			"(the name of the 'A' matrix)")
		stop(msg)
	}	
	if (typeof(S) != "character") {
		msg <- paste("argument 'S' is not a string",
			"(the name of the 'S' matrix)")
		stop(msg)
	}
	if (typeof(F) != "character") {
		msg <- paste("argument 'F' is not a string",
			"(the name of the 'F' matrix)")
		stop(msg)
	}
	if (!(single.na(M) || typeof(M) == "character")) {
		msg <- paste("argument M is not a string",
			"(the name of the 'M' matrix)")
		stop(msg)
	}
	if (is.na(M)) M <- as.integer(NA)
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
	threshnames <- checkThreshnames(threshnames)
	return(new("MxExpectationRAM", A, S, F, M, dimnames, thresholds, threshnames))
}

displayMxExpectationRAM <- function(expectation) {
	cat("MxExpectationRAM", omxQuotes(expectation@name), '\n')
	if (length(expectation@container)) {
		cat("$container :", omxQuotes(expectation@container), '\n')
	}
	if (length(expectation@submodels)) {
		cat("$submodels :", omxQuotes(expectation@submodels), '\n')
	}
	cat("$A :", omxQuotes(expectation@A), '\n')
	cat("$S :", omxQuotes(expectation@S), '\n')
	cat("$F :", omxQuotes(expectation@F), '\n')
	if (is.na(expectation@M)) {
		cat("$M :", expectation@M, '\n')
	} else {
		cat("$M :", omxQuotes(expectation@M), '\n')
	}
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
	invisible(expectation)
}

setMethod("print", "MxExpectationRAM", function(x,...) { 
	displayMxExpectationRAM(x) 
})

setMethod("show", "MxExpectationRAM", function(object) { 
	displayMxExpectationRAM(object) 
})
