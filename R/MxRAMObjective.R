#
#   Copyright 2007-2012 The OpenMx Project
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


setClass(Class = "MxRAMObjective",
	representation = representation(
		A = "MxCharOrNumber",
		S = "MxCharOrNumber",
		F = "MxCharOrNumber",
		M = "MxCharOrNumber",
		thresholds = "MxCharOrNumber",
		dims = "character",
		vector = "logical",
		definitionVars = "list",
		dataColumns = "numeric",
		thresholdColumns = "numeric",
		thresholdLevels = "numeric",
		depth = "integer",
		threshnames = "character",
		usePPML = "logical",
		ppmlData = "MxData"),
	contains = "MxBaseObjective")

setMethod("initialize", "MxRAMObjective",
	function(.Object, A, S, F, M, dims, thresholds, vector, threshnames,
		data = as.integer(NA), name = 'objective') {
		.Object@name <- name
		.Object@A <- A
		.Object@S <- S
		.Object@F <- F
		.Object@M <- M
		.Object@data <- data
		.Object@dims <- dims
		.Object@thresholds <- thresholds
		.Object@vector <- vector
		.Object@definitionVars <- list()
		.Object@threshnames <- threshnames
		.Object@usePPML <- FALSE
		return(.Object)
	}
)

setMethod("genericObjDependencies", signature("MxRAMObjective"),
	function(.Object, dependencies) {
	sources <- c(.Object@A, .Object@S, .Object@F, .Object@M, .Object@thresholds)
	sources <- sources[!is.na(sources)]
	dependencies <- imxAddDependency(sources, .Object@name, dependencies)
	return(dependencies)
})

setMethod("genericObjInitialMatrix", "MxRAMObjective",
	function(.Object, flatModel) {
		flatObjective <- flatModel@objectives[[.Object@name]]
		if (flatObjective@vector == FALSE) {
			return(matrix(as.double(NA), 1, 1))
		} else {
			modelname <- imxReverseIdentifier(flatModel, flatObjective@name)[[1]]
			name <- flatObjective@name
			if(is.na(flatObjective@data)) {
				msg <- paste("The RAM objective",
				"does not have a dataset associated with it in model",
				omxQuotes(modelname))
				stop(msg, call. = FALSE)
			}
			mxDataObject <- flatModel@datasets[[flatObjective@data]]
			if (mxDataObject@type != 'raw') {
				msg <- paste("The dataset associated with the RAM objective", 
					"in model", omxQuotes(modelname), "is not raw data.")
				stop(msg, call. = FALSE)
			}
			rows <- nrow(mxDataObject@observed)
			return(matrix(as.double(NA), rows, 1))
		}
})

setMethod("genericObjFunNamespace", signature("MxRAMObjective"), 
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

setMethod("genericObjRename", signature("MxRAMObjective"),
	function(.Object, oldname, newname) {
		.Object@A <- renameReference(.Object@A, oldname, newname)
		.Object@S <- renameReference(.Object@S, oldname, newname)
		.Object@F <- renameReference(.Object@F, oldname, newname)
		.Object@M <- renameReference(.Object@M, oldname, newname)
		.Object@data <- renameReference(.Object@data, oldname, newname)
		.Object@thresholds <- sapply(.Object@thresholds, renameReference, oldname, newname)		
		return(.Object)
})

setMethod("genericObjFunConvert", signature("MxRAMObjective", "MxFlatModel"), 
	function(.Object, flatModel, model, labelsData, defVars, dependencies) {
		modelname <- imxReverseIdentifier(model, .Object@name)[[1]]	
		name <- .Object@name
		aMatrix <- .Object@A
		sMatrix <- .Object@S
		fMatrix <- .Object@F
		mMatrix <- .Object@M
		data <- .Object@data
		if(is.na(data)) {
			msg <- paste("The RAM objective",
				"does not have a dataset associated with it in model",
				omxQuotes(modelname))
			stop(msg, call. = FALSE)
		}
		mxDataObject <- flatModel@datasets[[.Object@data]]
		if(!is.na(mMatrix) && single.na(mxDataObject@means) && mxDataObject@type != "raw") {
			msg <- paste("The RAM objective",
				"has an expected means vector but",
				"no observed means vector in model",
				omxQuotes(modelname))
			stop(msg, call. = FALSE)
		}
		if(!single.na(mxDataObject@means) && is.null(flatModel[[mMatrix]])) {
			msg <- paste("The RAM objective",
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
				"with the RAM objective function in model", 
				omxQuotes(modelname), "does not contain dimnames.")
				stop(msg, call. = FALSE)	
			}
			meanRows <- means[[1]]
			meanCols <- means[[2]]
			if (!is.null(meanRows) && length(meanRows) > 1) {
				msg <- paste("The M matrix associated",
				"with the RAM objective in model", 
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
		if (mxDataObject@type == 'raw') {
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

updateRAMdimnames <- function(flatObjective, job, flatJob, modelname) {
	fMatrixName <- flatObjective@F
	mMatrixName <- flatObjective@M
	if (is.na(mMatrixName)) {
		mMatrix <- NA
	} else {
		mMatrix <- job[[mMatrixName]]
	}
	fMatrix <- job[[fMatrixName]]
	if (is.null(fMatrix)) {
		stop(paste("Unknown F matrix name", 
			omxQuotes(simplifyName(fMatrixName, modelname)),
			"detected in the objective function",
			"of model", omxQuotes(modelname)), call. = FALSE)
	}
	dims <- flatObjective@dims
	if (!is.null(dimnames(fMatrix)) && !single.na(dims) && 
		!identical(dimnames(fMatrix)[[2]], dims)) {
		msg <- paste("The F matrix associated",
			"with the RAM objective in model", 
			omxQuotes(modelname), "contains dimnames and",
			"the objective function has specified dimnames")
		stop(msg, call.=FALSE)		
	}
	if (is.null(dimnames(fMatrix)) && !single.na(dims)) {
		fMatrixFlat <- flatJob[[fMatrixName]]
		dimnames(fMatrix) <- list(c(), dims)
		dimnames(fMatrixFlat) <- list(c(), dims)
		job[[fMatrixName]] <- fMatrix
		flatJob[[fMatrixName]] <- fMatrixFlat
	}
	if (!isS4(mMatrix) && (is.null(mMatrix) || is.na(mMatrix))) return(list(job, flatJob))
	if (!is.null(dimnames(mMatrix)) && !single.na(dims) &&
		!identical(dimnames(mMatrix), list(NULL, dims))) {
		msg <- paste("The M matrix associated",
			"with the RAM objective in model", 
			omxQuotes(modelname), "contains dimnames and",
			"the objective function has specified dimnames")
		stop(msg, call.=FALSE)	
	}
	if (is.null(dimnames(mMatrix)) && !single.na(dims)) {
		mMatrixFlat <- flatJob[[mMatrixName]]
		dimnames(mMatrix) <- list(NULL, dims)
		dimnames(mMatrixFlat) <- list(NULL, dims)
		job[[mMatrixName]] <- mMatrix
		flatJob[[mMatrixName]] <- mMatrixFlat
	}
	return(list(job, flatJob))
}

setMethod("genericObjModelConvert", "MxRAMObjective",
	function(.Object, job, model, namespace, labelsData, flatJob) {
		cache <- list()
		if(is.na(.Object@data)) {
			msg <- paste("The RAM objective",
				"does not have a dataset associated with it in model",
				omxQuotes(model@name))
			stop(msg, call.=FALSE)
		}
		
		tuple <- evaluateMxObject(.Object@A, flatJob, labelsData, cache)
		Amatrix <- tuple[[1]]
		cache <- tuple[[2]]
		tuple <- evaluateMxObject(.Object@S, flatJob, labelsData, cache)
		Smatrix <- tuple[[1]]
		cache <- tuple[[2]]
		if (!identical(dim(Amatrix), dim(Smatrix))) {
				msg <- paste("The RAM objective",
					"in model", omxQuotes(model@name), "has an A matrix",
					"with dimensions", nrow(Amatrix), 'x', ncol(Amatrix),
					"and a S matrix with dimensions", nrow(Smatrix), 'x',
					ncol(Smatrix))
				stop(msg, call.=FALSE)
		}

		# Simplest possible check for PPML applicability, imxTransformModelPPML
		# will check beyond that
		ppmlModelOption <- model@options$UsePPML
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

		pair <- updateRAMdimnames(.Object, job, flatJob, model@name)
		job <- pair[[1]]
		flatJob <- pair[[2]]
		if (flatJob@datasets[[.Object@data]]@type != 'raw') {
			if (.Object@vector) {
				msg <- paste("The RAM objective",
					"in model", omxQuotes(model@name), "has specified",
					"'vector' = TRUE, but the observed data is not raw data")
				stop(msg, call.=FALSE)
			}
			
			if (enablePPML) {
				oldName <- model@name
				model <- imxTransformModelPPML(model)
				# Make sure PPML was applied successfully
				if ( !is.null(model@options$UsePPML) ) {
					if ( model@options$UsePPML == "Split" ) {
						job[[oldName]] <- model
						
						job@.newobjects <- TRUE
						job@.newobjective <- TRUE
						job@.newtree <- TRUE
						return(list(job,flatJob))
					} else if ( model@options$UsePPML == "Solved" || model@options$UsePPML == "PartialSolved" ) { # TODO: PartialSolved
						job[[oldName]] <- model
						
						pair <- updateRAMdimnames(.Object, job, flatJob, model@name)
						job <- pair[[1]]
					
						job@.newobjects <- TRUE
						job@.newobjective <- FALSE
						job@.newtree <- FALSE
						return(list(job, flatJob))
					}
				} 
			}
			
			job@.newobjects <- FALSE
			job@.newobjective <- FALSE
			job@.newtree <- FALSE
			return(list(job, flatJob))
		}
		if (is.na(.Object@M) || is.null(job[[.Object@M]])) {
			msg <- paste("The RAM objective",
				"has raw data but is missing",
				"an expected means vector in model",
				omxQuotes(model@name))
			stop(msg, call.=FALSE)
		}
		
		if (enablePPML) {
			oldName <- model@name
			model <- imxTransformModelPPML(model)
			# Make sure PPML was applied successfully
			
			if ( !is.null(model@options$UsePPML) ) {
				if ( model@options$UsePPML == "Split" ) {
					job[[oldName]] <- model
					
					job@.newobjects <- TRUE
					job@.newobjective <- TRUE
					job@.newtree <- TRUE
					return(list(job,flatJob))
				} else if ( model@options$UsePPML == "Solved" || model@options$UsePPML == "PartialSolved" ) { # TODO: PartialSolved
					job[[oldName]] <- model
					
					pair <- updateRAMdimnames(.Object, job, flatJob, model@name)
					job <- pair[[1]]
					
					job@.newobjects <- TRUE
					job@.newobjective <- FALSE
					job@.newtree <- FALSE
					return(list(job, flatJob))
				}
			} 
		}
		
		pair <- updateThresholdDimnames(.Object, job, flatJob, model@name)
		job <- pair[[1]]
		flatJob <- pair[[2]]
		precision <- "Function precision"
		if(!single.na(.Object@thresholds) && 
			is.null(job@options[[precision]])) {
			job <- mxOption(job, precision, 1e-9)
			flatJob <- mxOption(flatJob, precision, 1e-9)
		}
		job@.newobjects <- FALSE
		job@.newobjective <- FALSE
		job@.newtree <- FALSE
		return(list(job, flatJob))
	}
)

imxSimpleRAMPredicate <- function(model) {
	if (is.null(model$objective) || !is(model$objective, "MxRAMObjective")) {
		return(FALSE)
	}
	nameA <- model$objective@A
	nameS <- model$objective@S
	A <- model[[nameA]]
	S <- model[[nameS]]
	if (is.null(A) || is.null(S)) {
		return(FALSE)
	}
	return(is(A, "MxMatrix") && is(S, "MxMatrix"))
}

mxRAMObjective <- function(A, S, F, M = NA, dimnames = NA, thresholds = NA, vector = FALSE,
									threshnames = dimnames) {
	if (missing(A) || typeof(A) != "character") {
		msg <- paste("argument 'A' is not a string",
			"(the name of the 'A' matrix)")
		stop(msg)
	}	
	if (missing(S) || typeof(S) != "character") {
		msg <- paste("argument 'S' is not a string",
			"(the name of the 'S' matrix)")
		stop(msg)
	}
	if (missing(F) || typeof(F) != "character") {
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
	if (single.na(threshnames)) threshnames <- as.character(NA)
	if (!is.vector(dimnames) || typeof(dimnames) != 'character') {
		stop("Dimnames argument is not a character vector")
	}
	if (!is.vector(threshnames) || typeof(threshnames) != 'character') {
		stop("'threshnames' argument is not a character vector")
	}
	if (length(thresholds) != 1) {
		stop("Thresholds argument must be a single matrix or algebra name")
	}
	if (length(dimnames) == 0) {
		stop("Dimnames argument cannot be an empty vector")
	}
	if (length(threshnames) == 0) {
		stop("'threshnames' argument cannot be an empty vector")
	}
	if (length(dimnames) > 1 && any(is.na(dimnames))) {
		stop("NA values are not allowed for dimnames vector")
	}
	if (length(threshnames) > 1 && any(is.na(threshnames))) {
		stop("NA values are not allowed for 'threshnames' vector")
	}
	if (length(vector) > 1 || typeof(vector) != "logical") {
		stop("Vector argument is not a logical value")
	}
	return(new("MxRAMObjective", A, S, F, M, dimnames, thresholds, vector, threshnames))
}

displayRAMObjective <- function(objective) {
	cat("MxRAMObjective", omxQuotes(objective@name), '\n')
	cat("@A :", omxQuotes(objective@A), '\n')
	cat("@S :", omxQuotes(objective@S), '\n')
	cat("@F :", omxQuotes(objective@F), '\n')
	if (is.na(objective@M)) {
		cat("@M :", objective@M, '\n')
	} else {
		cat("@M :", omxQuotes(objective@M), '\n')
	}
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
	cat("@info$likelihoods: ", length(objective@info$likelihoods) > 0, '\n')
	if (length(objective@result) == 0) {
		cat("@result: (not yet computed) ")
	} else {
		cat("@result:\n")
	}
	print(objective@result)
	if (length(objective@info$expCov) == 0) {
		cat("@info$expCov: ")
	} else {
		cat("@info$expCov:\n")
	}
	print(objective@info$expCov)
	if (length(objective@info$expMean) == 0) {
		cat("@info$expMean: ")
	} else {
		cat("@info$expMean:\n")
	}
	print(objective@info$expMean)
	invisible(objective)
}

setMethod("print", "MxRAMObjective", function(x,...) { 
	displayRAMObjective(x) 
})

setMethod("show", "MxRAMObjective", function(object) { 
	displayRAMObjective(object) 
})
