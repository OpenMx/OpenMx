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


# **DONE**
setClass(Class = "MxLISRELObjective",
	representation = representation(
		LX = "MxCharOrNumber",
		LY = "MxCharOrNumber",
		BE = "MxCharOrNumber",
		GA = "MxCharOrNumber",
		PH = "MxCharOrNumber",
		PS = "MxCharOrNumber",
		TD = "MxCharOrNumber",
		TE = "MxCharOrNumber",
		TH = "MxCharOrNumber",
		thresholds = "MxCharOrNumber",
		dims = "character",
		vector = "logical",
		expCov = "matrix",
		expMean = "matrix",
		definitionVars = "list",
		dataColumns = "numeric",
		thresholdColumns = "numeric",
		thresholdLevels = "numeric",
		depth = "integer"),
	contains = "MxBaseObjective")

# **DONE**
setMethod("initialize", "MxLISRELObjective",
	function(.Object, LX, LY, BE, GA, PH, PS, TD, TE, TH, dims, thresholds, vector,
		data = as.integer(NA), name = 'objective') {
		.Object@name <- name
		.Object@LX <- LX
		.Object@LY <- LY
		.Object@BE <- BE
		.Object@GA <- GA
		.Object@PH <- PH
		.Object@PS <- PS
		.Object@TD <- TD
		.Object@TE <- TE
		.Object@TH <- TH
		.Object@data <- data
		.Object@dims <- dims
		.Object@thresholds <- thresholds
		.Object@vector <- vector
		.Object@definitionVars <- list()
		return(.Object)
	}
)

# Not sure about the first line here
# Should it be MxLISRELObjective or MxBaseObjective?
setMethod("genericObjModelConvert", "MxLISRELObjective",
	function(.Object, job, model, namespace, flatJob) {
		if(is.na(.Object@data)) {
			msg <- paste("The LISREL objective",
				"does not have a dataset associated with it in model",
				omxQuotes(model@name))
			stop(msg, call.=FALSE)
		}
#		pair <- updateRAMdimnames(.Object, job, flatJob, model@name)
#		job <- pair[[1]]
#		flatJob <- pair[[2]]
#		if (flatJob@datasets[[.Object@data]]@type != 'raw') {
#			if (.Object@vector) {
#				msg <- paste("The RAM objective",
#					"in model", omxQuotes(model@name), "has specified",
#					"'vector' = TRUE, but the observed data is not raw data")
#				stop(msg, call.=FALSE)
#			}
#			job@.newobjects <- FALSE
#			job@.newobjective <- FALSE
#			job@.newtree <- FALSE
#			return(list(job, flatJob))
#		}
#		if (is.na(.Object@M) || is.null(job[[.Object@M]])) {
#			msg <- paste("The RAM objective",
#				"has raw data but is missing",
#				"an expected means vector in model",
#				omxQuotes(model@name))
#			stop(msg, call.=FALSE)
#		}
#		pair <- updateThresholdDimnames(.Object, job, flatJob, model@name)
#		job <- pair[[1]]
#		flatJob <- pair[[2]]
#		precision <- "Function precision"
#		if(!single.na(.Object@thresholds) && 
#			is.null(job@options[[precision]])) {
#			job <- mxOption(job, precision, 1e-9)
#			flatJob <- mxOption(flatJob, precision, 1e-9)
#		}
		job@.newobjects <- FALSE
		job@.newobjective <- FALSE
		job@.newtree <- FALSE
		return(list(job, flatJob))
	}
)

# **DONE**
setMethod("genericObjFunNamespace", signature("MxLISRELObjective"), 
	function(.Object, modelname, namespace) {
		.Object@name <- imxIdentifier(modelname, .Object@name)
		.Object@LX <- imxConvertIdentifier(.Object@LX, modelname, namespace)
		.Object@LY <- imxConvertIdentifier(.Object@LY, modelname, namespace)
		.Object@BE <- imxConvertIdentifier(.Object@BE, modelname, namespace)
		.Object@GA <- imxConvertIdentifier(.Object@GA, modelname, namespace)
		.Object@PH <- imxConvertIdentifier(.Object@PH, modelname, namespace)
		.Object@PS <- imxConvertIdentifier(.Object@PS, modelname, namespace)
		.Object@TD <- imxConvertIdentifier(.Object@TD, modelname, namespace)
		.Object@TE <- imxConvertIdentifier(.Object@TE, modelname, namespace)
		.Object@TH <- imxConvertIdentifier(.Object@TH, modelname, namespace)
		.Object@data <- imxConvertIdentifier(.Object@data, modelname, namespace)
		.Object@thresholds <- sapply(.Object@thresholds, imxConvertIdentifier, modelname, namespace)
		return(.Object)
	}
)


# **DONE**
setMethod("genericObjFunConvert", signature("MxLISRELObjective", "MxFlatModel"), 
	function(.Object, flatModel, model, defVars) {
		modelname <- imxReverseIdentifier(model, .Object@name)[[1]]	
		name <- .Object@name
		lxMatrix <- .Object@LX
		lyMatrix <- .Object@LY
		beMatrix <- .Object@BE
		gaMatrix <- .Object@GA
		phMatrix <- .Object@PH
		psMatrix <- .Object@PS
		tdMatrix <- .Object@TD
		teMatrix <- .Object@TE
		thMatrix <- .Object@TH
		data <- .Object@data
		if(is.na(data)) {
			msg <- paste("The LISREL objective",
				"does not have a dataset associated with it in model",
				omxQuotes(modelname))
			stop(msg, call. = FALSE)
		}
		mxDataObject <- flatModel@datasets[[.Object@data]]
#		if(!is.na(mMatrix) && single.na(mxDataObject@means) && mxDataObject@type != "raw") {
#			msg <- paste("The RAM objective",
#				"has an expected means vector but",
#				"no observed means vector in model",
#				omxQuotes(modelname))
#			stop(msg, call. = FALSE)
#		}
#		if(!single.na(mxDataObject@means) && is.null(flatModel[[mMatrix]])) {
#			msg <- paste("The RAM objective",
#				"has an observed means vector but",
#				"no expected means vector in model",
#				omxQuotes(modelname))
#			stop(msg, call. = FALSE)		
#		}
		checkNumericData(mxDataObject)
		.Object@LX <- imxLocateIndex(flatModel, lxMatrix, name)
		.Object@LY <- imxLocateIndex(flatModel, lyMatrix, name)
		.Object@BE <- imxLocateIndex(flatModel, beMatrix, name)
		.Object@GA <- imxLocateIndex(flatModel, gaMatrix, name)
		.Object@PH <- imxLocateIndex(flatModel, phMatrix, name)
		.Object@PS <- imxLocateIndex(flatModel, psMatrix, name)
		.Object@TD <- imxLocateIndex(flatModel, tdMatrix, name)
		.Object@TE <- imxLocateIndex(flatModel, teMatrix, name)
		.Object@TH <- imxLocateIndex(flatModel, thMatrix, name)
		.Object@data <- as.integer(imxLocateIndex(flatModel, data, name))
#		verifyObservedNames(mxDataObject@observed, mxDataObject@means, mxDataObject@type, flatModel, modelname, "RAM")
#		fMatrix <- flatModel[[fMatrix]]@values
#		if (is.null(dimnames(fMatrix))) {
#			msg <- paste("The F matrix of model",
#				omxQuotes(modelname), "does not contain dimnames")
#			stop(msg, call. = FALSE)
#		}
#		if (is.null(dimnames(fMatrix)[[2]])) {
#			msg <- paste("The F matrix of model",
#				omxQuotes(modelname), "does not contain colnames")
#			stop(msg, call. = FALSE)
#		}
#		mMatrix <- flatModel[[mMatrix]]		
#		if(!is.null(mMatrix)) {
#			means <- dimnames(mMatrix)
#			if (is.null(means)) {
#				msg <- paste("The M matrix associated",
#				"with the RAM objective function in model", 
#				omxQuotes(modelname), "does not contain dimnames.")
#				stop(msg, call. = FALSE)	
#			}
#			meanRows <- means[[1]]
#			meanCols <- means[[2]]
#			if (!is.null(meanRows) && length(meanRows) > 1) {
#				msg <- paste("The M matrix associated",
#				"with the RAM objective in model", 
#				omxQuotes(modelname), "is not a 1 x N matrix.")
#				stop(msg, call. = FALSE)
#			}
#			if (!identical(dimnames(fMatrix)[[2]], meanCols)) {
#				msg <- paste("The column names of the F matrix",
#					"and the column names of the M matrix",
#					"in model", 
#					omxQuotes(modelname), "do not contain identical",
#					"names.")
#				stop(msg, call. = FALSE)
#			}
#		}
#		translatedNames <- fMatrixTranslateNames(fMatrix, modelname)
#		.Object@depth <- generateRAMDepth(flatModel, aMatrix, model@options)
#		if (mxDataObject@type == 'raw') {
#			threshName <- .Object@thresholds
#			checkNumberOrdinalColumns(mxDataObject)
#			.Object@definitionVars <- imxFilterDefinitionVariables(defVars, data)
#			.Object@dataColumns <- generateDataColumns(flatModel, translatedNames, data)
#			verifyThresholds(flatModel, model, data, translatedNames, threshName)
#			.Object@thresholds <- imxLocateIndex(flatModel, threshName, name)
#			retval <- generateThresholdColumns(flatModel, model, translatedNames, data, threshName)
#			.Object@thresholdColumns <- retval[[1]]
#			.Object@thresholdLevels <- retval[[2]]
#			if (length(mxDataObject@observed) == 0) {
#				.Object@data <- as.integer(NA)
#			}
#			if (single.na(.Object@dims)) {
#				.Object@dims <- translatedNames
#			}
#		} else {
#			if (!identical(translatedNames, rownames(mxDataObject@observed))) {
#				msg <- paste("The names of the manifest",
#					"variables in the F matrix of model",
#					omxQuotes(modelname), "does not match the",
#					"dimnames of the observed covariance matrix")
#				stop(msg, call. = FALSE)
#			}
#		}
		return(.Object)
	}
)


# **DONE**
setMethod("genericObjDependencies", signature("MxLISRELObjective"),
	function(.Object, dependencies) {
	sources <- c(.Object@LX, .Object@LY, .Object@BE, .Object@GA, .Object@PH, .Object@PS, .Object@TD, .Object@TE, .Object@TH, .Object@thresholds)
	sources <- sources[!is.na(sources)]
	dependencies <- imxAddDependency(sources, .Object@name, dependencies)
	return(dependencies)
	}
)


# **DONE**
setMethod("genericObjRename", signature("MxLISRELObjective"),
	function(.Object, oldname, newname) {
		.Object@LX <- renameReference(.Object@LX, oldname, newname)
		.Object@LY <- renameReference(.Object@LY, oldname, newname)
		.Object@BE <- renameReference(.Object@BE, oldname, newname)
		.Object@GA <- renameReference(.Object@GA, oldname, newname)
		.Object@PH <- renameReference(.Object@PH, oldname, newname)
		.Object@PS <- renameReference(.Object@PS, oldname, newname)
		.Object@TD <- renameReference(.Object@TD, oldname, newname)
		.Object@TE <- renameReference(.Object@TE, oldname, newname)
		.Object@TH <- renameReference(.Object@TH, oldname, newname)
		.Object@data <- renameReference(.Object@data, oldname, newname)
		.Object@thresholds <- sapply(.Object@thresholds, renameReference, oldname, newname)		
		return(.Object)
	}
)


# **DONE**
mxLISRELObjective <- function(LX, LY, BE, GA, PH, PS, TD, TE, TH, dimnames = NA, thresholds = NA, vector = FALSE) {
	if (missing(LX) || typeof(LX) != "character") {
		msg <- paste("argument 'LX' is not a string",
			"(the name of the 'LX' matrix)")
		stop(msg)
	}
	if (missing(LY) || typeof(LY) != "character") {
		msg <- paste("argument 'LY' is not a string",
			"(the name of the 'LY' matrix)")
		stop(msg)
	}
	if (missing(BE) || typeof(BE) != "character") {
		msg <- paste("argument 'BE' is not a string",
			"(the name of the 'BE' matrix)")
		stop(msg)
	}
	if (missing(GA) || typeof(GA) != "character") {
		msg <- paste("argument 'GA' is not a string",
			"(the name of the 'GA' matrix)")
		stop(msg)
	}
	if (missing(PH) || typeof(PH) != "character") {
		msg <- paste("argument 'PH' is not a string",
			"(the name of the 'PH' matrix)")
		stop(msg)
	}
	if (missing(PS) || typeof(PS) != "character") {
		msg <- paste("argument 'PS' is not a string",
			"(the name of the 'PS' matrix)")
		stop(msg)
	}
	if (missing(TD) || typeof(TD) != "character") {
		msg <- paste("argument 'TD' is not a string",
			"(the name of the 'TD' matrix)")
		stop(msg)
	}
	if (missing(TE) || typeof(TE) != "character") {
		msg <- paste("argument 'TE' is not a string",
			"(the name of the 'TE' matrix)")
		stop(msg)
	}
	if (missing(TH) || typeof(TH) != "character") {
		msg <- paste("argument 'TH' is not a string",
			"(the name of the 'TH' matrix)")
		stop(msg)
	}
#	if (!(single.na(M) || typeof(M) == "character")) {
#		msg <- paste("argument M is not a string",
#			"(the name of the 'M' matrix)")
#		stop(msg)
#	}
#	if (is.na(M)) M <- as.integer(NA)
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
	if (length(vector) > 1 || typeof(vector) != "logical") {
		stop("Vector argument is not a logical value")
	}
	return(new("MxLISRELObjective", LX, LY, BE, GA, PH, PS, TD, TE, TH, dimnames, thresholds, vector))
}


# **DONE**
displayLISRELObjective <- function(objective) {
	cat("MxLISRELObjective", omxQuotes(objective@name), '\n')
	cat("@LX :", omxQuotes(objective@LX), '\n')
	cat("@LY :", omxQuotes(objective@LY), '\n')
	cat("@BE :", omxQuotes(objective@BE), '\n')
	cat("@GA :", omxQuotes(objective@GA), '\n')
	cat("@PH :", omxQuotes(objective@PH), '\n')
	cat("@PS :", omxQuotes(objective@PS), '\n')
	cat("@TD :", omxQuotes(objective@TD), '\n')
	cat("@TE :", omxQuotes(objective@TE), '\n')
	cat("@TH :", omxQuotes(objective@TH), '\n')
#	if (is.na(objective@M)) {
#		cat("@M :", objective@M, '\n')
#	} else {
#		cat("@M :", omxQuotes(objective@M), '\n')
#	}
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
	if (length(objective@result) == 0) {
		cat("@result: (not yet computed) ")
	} else {
		cat("@result:\n")
	}
	print(objective@result)
	if (length(objective@expCov) == 0) {
		cat("@expCov: ")
	} else {
		cat("@expCov:\n")
	}
	print(objective@expCov)
	if (length(objective@expMean) == 0) {
		cat("@expMean: ")
	} else {
		cat("@expMean:\n")
	}
	print(objective@expMean)
	invisible(objective)
}


# **DONE**
setMethod("print", "MxLISRELObjective", function(x,...) { 
	displayLISRELObjective(x) 
})


# **DONE**
setMethod("show", "MxLISRELObjective", function(object) { 
	displayLISRELObjective(object) 
})





#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
# BEGIN SECTION OF THINGS I DO NOT THINK I NEED



#setMethod("genericObjInitialMatrix", "MxRAMObjective",
#	function(.Object, flatModel) {
#		flatObjective <- flatModel@objectives[[.Object@name]]
#		if (flatObjective@vector == FALSE) {
#			return(matrix(as.double(NA), 1, 1))
#		} else {
#			modelname <- imxReverseIdentifier(flatModel, flatObjective@name)[[1]]
#			name <- flatObjective@name
#			if(is.na(flatObjective@data)) {
#				msg <- paste("The RAM objective",
#				"does not have a dataset associated with it in model",
#				omxQuotes(modelname))
#				stop(msg, call. = FALSE)
#			}
#			mxDataObject <- flatModel@datasets[[flatObjective@data]]
#			if (mxDataObject@type != 'raw') {
#				msg <- paste("The dataset associated with the RAM objective", 
#					"in model", omxQuotes(modelname), "is not raw data.")
#				stop(msg, call. = FALSE)
#			}
#			rows <- nrow(mxDataObject@observed)
#			return(matrix(as.double(NA), rows, 1))
#		}
#})
#
#
#
#setMethod("genericObjReadAttributes", signature("MxRAMObjective"),
#	function(.Object, values) {
#		.Object@expCov <- attr(values, "expCov", exact = TRUE)
#		.Object@expMean <- attr(values, "expMean", exact = TRUE)
#		dimnames(values) <- dimnames(.Object)
#		attr(values, "expCov") <- NULL
#		attr(values, "expMean") <- NULL
#		.Object@result <- values
#		return(.Object)
#})
#
#generateRAMDepth <- function(flatModel, aMatrixName, modeloptions) {
#	mxObject <- flatModel[[aMatrixName]]
#	if (!is(mxObject, "MxMatrix")) {
#		return(as.integer(NA))
#	}
#	if (identical(modeloptions[['RAM Inverse Optimization']], "No")) {
#		return(as.integer(NA))
#	}
#	if (is.null(modeloptions[['RAM Inverse Optimization']]) &&
#		identical(getOption('mxOptions')[['RAM Inverse Optimization']], "No")) {
#		return(as.integer(NA))
#	}	
#	maxdepth <- modeloptions[['RAM Max Depth']]
#	if (is.null(maxdepth) || (length(maxdepth) != 1) ||
#		is.na(maxdepth) || !is.numeric(maxdepth) || maxdepth < 0) {
#		maxdepth <- nrow(mxObject) - 1
#	}
#	return(omxGetRAMDepth(mxObject, maxdepth))
#}
#
#omxGetRAMDepth <- function(A, maxdepth = nrow(A) - 1) {
#	mxObject <- A
#	aValues <- matrix(0, nrow(mxObject), ncol(mxObject))
#	defvars <- apply(mxObject@labels, c(1,2), imxIsDefinitionVariable)
#	squarebrackets <- apply(mxObject@labels, c(1,2), hasSquareBrackets)
#	aValues[mxObject@free] <- 1
#	aValues[mxObject@values != 0] <- 1
#	aValues[defvars] <- 1
#	aValues[squarebrackets] <- 1
#	return(generateDepthHelper(aValues, aValues, 0, maxdepth))
#}
#
#generateDepthHelper <- function(aValues, currentProduct, depth, maxdepth) {
#	if (depth > maxdepth) {
#		return(as.integer(NA))
#	}
#	if (all(currentProduct == 0)) { 
#		return(as.integer(depth))
#	} else {
#		return(generateDepthHelper(aValues, currentProduct %*% aValues, depth + 1, maxdepth))
#	}
#}
#
#fMatrixTranslateNames <- function(fMatrix, modelName) {
#	retval <- character()
#	colNames <- dimnames(fMatrix)[[2]]
#	for(i in 1:nrow(fMatrix)) {
#		irow <- fMatrix[i,]
#		matches <- which(irow == 1)
#		if (length(matches) != 1) {
#			err <- paste("The model",
#				omxQuotes(modelName), "does not contain",
#				"a valid F matrix")
#			stop(err, call. = FALSE)
#		}
#		retval[[i]] <- colNames[[matches[[1]]]]
#	}
#	return(retval)
#}
#
#updateRAMdimnames <- function(flatObjective, job, flatJob, modelname) {
#	fMatrixName <- flatObjective@F
#	mMatrixName <- flatObjective@M
#	if (is.na(mMatrixName)) {
#		mMatrix <- NA
#	} else {
#		mMatrix <- job[[mMatrixName]]
#	}
#	fMatrix <- job[[fMatrixName]]
#	if (is.null(fMatrix)) {
#		stop(paste("Unknown F matrix name", 
#			omxQuotes(simplifyName(fMatrixName, modelname)),
#			"detected in the objective function",
#			"of model", omxQuotes(modelname)), call. = FALSE)
#	}
#	dims <- flatObjective@dims
#	if (!is.null(dimnames(fMatrix)) && !single.na(dims) && 
#		!identical(dimnames(fMatrix)[[2]], dims)) {
#		msg <- paste("The F matrix associated",
#			"with the RAM objective in model", 
#			omxQuotes(modelname), "contains dimnames and",
#			"the objective function has specified dimnames")
#		stop(msg, call.=FALSE)		
#	}
#	if (is.null(dimnames(fMatrix)) && !single.na(dims)) {
#		fMatrixFlat <- flatJob[[fMatrixName]]
#		dimnames(fMatrix) <- list(c(), dims)
#		dimnames(fMatrixFlat) <- list(c(), dims)
#		job[[fMatrixName]] <- fMatrix
#		flatJob[[fMatrixName]] <- fMatrixFlat
#	}
#	if (!isS4(mMatrix) && (is.null(mMatrix) || is.na(mMatrix))) return(list(job, flatJob))
#	if (!is.null(dimnames(mMatrix)) && !single.na(dims) &&
#		!identical(dimnames(mMatrix), list(NULL, dims))) {
#		msg <- paste("The M matrix associated",
#			"with the RAM objective in model", 
#			omxQuotes(modelname), "contains dimnames and",
#			"the objective function has specified dimnames")
#		stop(msg, call.=FALSE)	
#	}
#	if (is.null(dimnames(mMatrix)) && !single.na(dims)) {
#		mMatrixFlat <- flatJob[[mMatrixName]]
#		dimnames(mMatrix) <- list(NULL, dims)
#		dimnames(mMatrixFlat) <- list(NULL, dims)
#		job[[mMatrixName]] <- mMatrix
#		flatJob[[mMatrixName]] <- mMatrixFlat
#	}
#	return(list(job, flatJob))
#}
# END SECTION OF THINGS I DO NO THINK I NEED
