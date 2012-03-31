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

#--------------------------------------------------------------------
# Author: Michael D. Hunter
# Filename: MxLISRELObjective.R
#--------------------------------------------------------------------

#--------------------------------------------------------------------
# Revision History
#   Mon Feb 20 13:03:21 Central Standard Time 2012 -- Michael Hunter added means


#--------------------------------------------------------------------
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
		TX = "MxCharOrNumber",
		TY = "MxCharOrNumber",
		KA = "MxCharOrNumber",
		AL = "MxCharOrNumber",
		thresholds = "MxCharOrNumber",
		dims = "character",
		vector = "logical",
		expCov = "matrix",
		expMean = "matrix",
		definitionVars = "list",
		dataColumns = "numeric", #Used in FIML to grab the correct data elements
		thresholdColumns = "numeric", #Used in FIML
		thresholdLevels = "numeric", # Used in FIML
		threshnames = "character",
		depth = "integer"), #Used to speed up I-A inverse in RAM, could be used to speed up I-B inverse in LISREL
	contains = "MxBaseObjective")

# **DONE**
setMethod("initialize", "MxLISRELObjective",
	function(.Object, LX, LY, BE, GA, PH, PS, TD, TE, TH, TX, TY, KA, AL, dims, thresholds, vector, threshnames,
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
		.Object@TX <- TX
		.Object@TY <- TY
		.Object@KA <- KA
		.Object@AL <- AL
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
	function(.Object, job, model, namespace, labelsData, flatJob) {
		if(is.na(.Object@data)) {
			msg <- paste("The LISREL objective",
				"does not have a dataset associated with it in model",
				omxQuotes(model@name))
			stop(msg, call.=FALSE)
		}
#		The code below is out of date.  See current MxRAMObjective for up to date info.
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
		.Object@TX <- imxConvertIdentifier(.Object@TX, modelname, namespace)
		.Object@TY <- imxConvertIdentifier(.Object@TY, modelname, namespace)
		.Object@KA <- imxConvertIdentifier(.Object@KA, modelname, namespace)
		.Object@AL <- imxConvertIdentifier(.Object@AL, modelname, namespace)
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
		txMatrix <- .Object@TX
		tyMatrix <- .Object@TY
		kaMatrix <- .Object@KA
		alMatrix <- .Object@AL
		data <- .Object@data
		if(is.na(data)) {
			msg <- paste("The LISREL objective",
				"does not have a dataset associated with it in model",
				omxQuotes(modelname))
			stop(msg, call. = FALSE)
		}
		mxDataObject <- flatModel@datasets[[.Object@data]]
		if(!is.na(txMatrix) && single.na(mxDataObject@means) && mxDataObject@type != "raw") {
			msg <- paste("The LISREL objective",
				"has an expected means vector, TX, but",
				"no observed means vector in model",
				omxQuotes(modelname))
			stop(msg, call. = FALSE)
		}
		if(!is.na(tyMatrix) && single.na(mxDataObject@means) && mxDataObject@type != "raw") {
			msg <- paste("The LISREL objective",
				"has an expected means vector, TY, but",
				"no observed means vector in model",
				omxQuotes(modelname))
			stop(msg, call. = FALSE)
		}
		if(!is.na(kaMatrix) && single.na(mxDataObject@means) && mxDataObject@type != "raw") {
			msg <- paste("The LISREL objective",
				"has an expected means vector, KA, but",
				"no observed means vector in model",
				omxQuotes(modelname))
			stop(msg, call. = FALSE)
		}
		if(!is.na(alMatrix) && single.na(mxDataObject@means) && mxDataObject@type != "raw") {
			msg <- paste("The LISREL objective",
				"has an expected means vector, AL, but",
				"no observed means vector in model",
				omxQuotes(modelname))
			stop(msg, call. = FALSE)
		}
		if(!single.na(mxDataObject@means) && (is.null(flatModel[[txMatrix]]) || is.null(flatModel[[tyMatrix]]) || is.null(flatModel[[kaMatrix]]) || is.null(flatModel[[alMatrix]]) )) {
			msg <- paste("The LISREL objective",
				"has an observed means vector but",
				"is missing some expected means vectors in model",
				omxQuotes(modelname))
			stop(msg, call. = FALSE)
		}
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
		.Object@TX <- imxLocateIndex(flatModel, txMatrix, name)
		.Object@TY <- imxLocateIndex(flatModel, tyMatrix, name)
		.Object@KA <- imxLocateIndex(flatModel, kaMatrix, name)
		.Object@AL <- imxLocateIndex(flatModel, alMatrix, name)
		.Object@data <- as.integer(imxLocateIndex(flatModel, data, name))
		
		# 
		verifyObservedNames(mxDataObject@observed, mxDataObject@means, mxDataObject@type, flatModel, modelname, "LISREL")
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
		txMatrix <- flatModel[[txMatrix]]
		tyMatrix <- flatModel[[tyMatrix]]
		kaMatrix <- flatModel[[kaMatrix]]
		alMatrix <- flatModel[[alMatrix]]
		lxMatrix <- flatModel[[lxMatrix]]
		lyMatrix <- flatModel[[lyMatrix]]
		if(!is.null(txMatrix) || !is.null(tyMatrix) || !is.null(kaMatrix) || !is.null(alMatrix) ) {
#			means <- dimnames(mMatrix)
#			if (is.null(means)) { #Check if the means have dimnames
#				msg <- paste("The M matrix associated",
#				"with the RAM objective function in model", 
#				omxQuotes(modelname), "does not contain dimnames.")
#				stop(msg, call. = FALSE)	
#			}
#			meanRows <- means[[1]]
#			meanCols <- means[[2]]
#			if (!is.null(meanRows) && length(meanRows) > 1) { #Check if means are a row/column vector
#				msg <- paste("The M matrix associated",
#				"with the RAM objective in model", 
#				omxQuotes(modelname), "is not a 1 x N matrix.")
#				stop(msg, call. = FALSE)
#			}
#			if (!identical(dimnames(fMatrix)[[2]], meanCols)) { #Check if means exactly match F matrix (including order)
#				msg <- paste("The column names of the F matrix",
#					"and the column names of the M matrix",
#					"in model", 
#					omxQuotes(modelname), "do not contain identical",
#					"names.")
#				stop(msg, call. = FALSE)
#			}
		}
		translatedNames <- c(dimnames(lyMatrix)[[1]], dimnames(lxMatrix)[[1]]) #fMatrixTranslateNames(fMatrix, modelname) #Rearrange the rownames of F to match the order of the columns
#		.Object@depth <- generateRAMDepth(flatModel, aMatrix, model@options)
		if (mxDataObject@type == 'raw') {
			threshName <- .Object@thresholds
			checkNumberOrdinalColumns(mxDataObject)
			.Object@definitionVars <- imxFilterDefinitionVariables(defVars, data)
			.Object@dataColumns <- generateDataColumns(flatModel, translatedNames, data)
			verifyThresholds(flatModel, model, data, translatedNames, threshName)
			.Object@thresholds <- imxLocateIndex(flatModel, threshName, name)
			retval <- generateThresholdColumns(flatModel, model, translatedNames, data, threshName)
			.Object@thresholdColumns <- retval[[1]]
			.Object@thresholdLevels <- retval[[2]]
			if (length(mxDataObject@observed) == 0) {
				.Object@data <- as.integer(NA)
			}
			if (single.na(.Object@dims)) {
				.Object@dims <- translatedNames
			}
		}# else {
#			if (!identical(translatedNames, rownames(mxDataObject@observed))) { #Check the rows of F match the obs Cov
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
	sources <- c(.Object@LX, .Object@LY, .Object@BE, .Object@GA, .Object@PH, .Object@PS, .Object@TD, .Object@TE, .Object@TH, .Object@TX, .Object@TY, .Object@KA, .Object@AL, .Object@thresholds)
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
		.Object@TX <- renameReference(.Object@TX, oldname, newname)
		.Object@TY <- renameReference(.Object@TY, oldname, newname)
		.Object@KA <- renameReference(.Object@KA, oldname, newname)
		.Object@AL <- renameReference(.Object@AL, oldname, newname)
		.Object@data <- renameReference(.Object@data, oldname, newname)
		.Object@thresholds <- sapply(.Object@thresholds, renameReference, oldname, newname)		
		return(.Object)
	}
)


# **DONE**
mxLISRELObjective <- function(LX, LY, BE, GA, PH, PS, TD, TE, TH, TX = NA, TY = NA, KA = NA, AL = NA, dimnames = NA, thresholds = NA, vector = FALSE, threshnames = dimnames) {
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
	if (!(single.na(TX) || typeof(TX) == "character")) {
		msg <- paste("argument TX is not a string",
			"(the name of the 'TX' matrix)")
		stop(msg)
	}
	if (is.na(TX)) TX <- as.integer(NA)
	if (!(single.na(TY) || typeof(TY) == "character")) {
		msg <- paste("argument TY is not a string",
			"(the name of the 'TY' matrix)")
		stop(msg)
	}
	if (is.na(TY)) TY <- as.integer(NA)
	if (!(single.na(KA) || typeof(KA) == "character")) {
		msg <- paste("argument KA is not a string",
			"(the name of the 'KA' matrix)")
		stop(msg)
	}
	if (is.na(KA)) KA <- as.integer(NA)
	if (!(single.na(AL) || typeof(AL) == "character")) {
		msg <- paste("argument AL is not a string",
			"(the name of the 'AL' matrix)")
		stop(msg)
	}
	if (is.na(AL)) AL <- as.integer(NA)
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
	return(new("MxLISRELObjective", LX, LY, BE, GA, PH, PS, TD, TE, TH, TX, TY, KA, AL, dimnames, thresholds, vector, threshnames))
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
	if (is.na(objective@TX)) {
		cat("@TX :", objective@TX, '\n')
	} else {
		cat("@TX :", omxQuotes(objective@TX), '\n')
	}
	if (is.na(objective@TY)) {
		cat("@TY :", objective@TY, '\n')
	} else {
		cat("@TY :", omxQuotes(objective@TY), '\n')
	}
	if (is.na(objective@KA)) {
		cat("@KA :", objective@KA, '\n')
	} else {
		cat("@KA :", omxQuotes(objective@KA), '\n')
	}
	if (is.na(objective@AL)) {
		cat("@AL :", objective@AL, '\n')
	} else {
		cat("@AL :", omxQuotes(objective@AL), '\n')
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
