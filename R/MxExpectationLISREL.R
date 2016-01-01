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

#--------------------------------------------------------------------
# Author: Michael D. Hunter
# Filename: MxLISRELObjective.R
#--------------------------------------------------------------------

#--------------------------------------------------------------------
# Revision History
#   Mon Feb 20 13:03:21 Central Standard Time 2012 -- Michael Hunter added means
#   Sat Apr 07 19:48:33 Central Daylight Time 2012 -- Michael Hunter added lots of error checking
# 


#--------------------------------------------------------------------
# **DONE**
setClass(Class = "MxExpectationLISREL",
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
	        numStats = "numeric",
		thresholds = "MxCharOrNumber",
		dims = "character",
		dataColumns = "numeric", #Used in FIML to grab the correct data elements
		thresholdColumns = "numeric", #Used in FIML
		thresholdLevels = "numeric", # Used in FIML
		threshnames = "character",
		depth = "integer"), #Used to speed up I-A inverse in RAM, could be used to speed up I-B inverse in LISREL
	contains = "BaseExpectationNormal")


#--------------------------------------------------------------------
# **DONE**
setMethod("initialize", "MxExpectationLISREL",
	function(.Object, LX, LY, BE, GA, PH, PS, TD, TE, TH, TX, TY, KA, AL, dims, thresholds, threshnames,
		data = as.integer(NA), name = 'expectation') {
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
		return(.Object)
	}
)


#--------------------------------------------------------------------
setMethod("genericExpConvertEntities", "MxExpectationLISREL",
	function(.Object, flatModel, namespace, labelsData) {
		if(is.na(.Object@data)) {
			modelname <- getModelName(.Object)
			msg <- paste("The LISREL expectation function",
				"does not have a dataset associated with it in model",
				omxQuotes(modelname))
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
		return(flatModel)
	}
)


#--------------------------------------------------------------------
# **DONE**
setMethod("qualifyNames", signature("MxExpectationLISREL"), 
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


#--------------------------------------------------------------------
# Helper functions used in genericObjFunConvert method

checkLISRELmeansHelper <- function(Lam, Mean, Latent, matrixname, lamname, modelname){
	if(Latent){
		varType <- "latent" #used in error messages
		checkInd <- 2  #used to check row or col names match (1 for rows, 2 for cols)
		checkStg <- "col"
	}
	else{
		varType <- "manifest"
		checkInd <- 1
		checkStg <- "row"
	}
	# Check that the means are non-null
	if(is.null(Mean)){
		msg <- paste("The LISREL expectation function",
			"has an observed means vector but",
			"is missing expected means vector",
			"for", varType, "variables",
			matrixname,  "in model",
			omxQuotes(modelname))
		stop(msg, call. = FALSE)
	}
	# Check that the means have dimnames
	meanDimnames <- dimnames(Mean)
	if(is.null(meanDimnames)){
		msg <- paste("The", matrixname, "matrix associated",
		"with the LISREL expectation function in model", 
		omxQuotes(modelname), "does not contain dimnames.")
		stop(msg, call. = FALSE)	
	}
	#Check if means are a column vector
	meanRownames <- meanDimnames[[1]]
	meanColnames <- meanDimnames[[2]]
	if (!is.null(meanColnames) && length(meanColnames) > 1) {
		msg <- paste("The", matrixname, "matrix associated",
		"with the LISREL expectation function in model", 
		omxQuotes(modelname), "is not an N x 1 matrix.")
		stop(msg, call. = FALSE)
	}
	#Check if means exactly match Lambda matrix (including order)
	if (!identical(dimnames(Lam)[[checkInd]], meanRownames)) {
		msg <- paste("The", checkStg, "names of the", lamname, "matrix",
			"and the row names of the", matrixname, "matrix",
			"in model", 
			omxQuotes(modelname), "do not contain identical",
			"names.")
		stop(msg, call. = FALSE)
	}
}


checkLISRELmeans <- function(Lam, ManMean, LatMean, X, modelname){
	if(X){
		manMeanMat <- 'TX'
		latMeanMat <- 'KA'
		lamMat <- 'LX'
	} else{
		manMeanMat <- 'TY'
		latMeanMat <- 'AL'
		lamMat <- 'LY'
	}
	checkLISRELmeansHelper(
		Lam=Lam,
		Mean=ManMean,
		Latent=FALSE,
		matrixname=manMeanMat,
		lamname= lamMat,
		modelname=modelname
	)
	checkLISRELmeansHelper(
		Lam=Lam,
		Mean=LatMean,
		Latent=TRUE,
		matrixname=latMeanMat,
		lamname= lamMat,
		modelname=modelname
	)
}


#--------------------------------------------------------------------
# **DONE**
# Note: Lots of error checking is done in this method
setMethod("genericExpFunConvert", signature("MxExpectationLISREL"), 
	function(.Object, flatModel, model, labelsData, dependencies) {
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
		beMatrix2 <- beMatrix #This is a placeholder for use with the I-BE inverse speedup
		# Check if the model has data
		if(is.na(data)) {
			msg <- paste("The LISREL expectation function",
				"does not have a dataset associated with it in model",
				omxQuotes(modelname))
			stop(msg, call. = FALSE)
		}
		#
		# if any of the names of TX, TY, KA, AL are not missing
		#	 then the model must have observed means or raw data
		mxDataObject <- flatModel@datasets[[.Object@data]]
		if(single.na(mxDataObject@means) && mxDataObject@type != "raw") {
			if(!is.na(txMatrix)) {
				msg <- paste("The LISREL expectation function",
					"has an expected means vector, TX, but",
					"no observed means vector in model",
					omxQuotes(modelname))
				stop(msg, call. = FALSE)
			}
			if(!is.na(tyMatrix)) {
				msg <- paste("The LISREL expectation function",
					"has an expected means vector, TY, but",
					"no observed means vector in model",
					omxQuotes(modelname))
				stop(msg, call. = FALSE)
			}
			if(!is.na(kaMatrix)) {
				msg <- paste("The LISREL expectation function",
					"has an expected means vector, KA, but",
					"no observed means vector in model",
					omxQuotes(modelname))
				stop(msg, call. = FALSE)
			}
			if(!is.na(alMatrix)) {
				msg <- paste("The LISREL expectation function",
					"has an expected means vector, AL, but",
					"no observed means vector in model",
					omxQuotes(modelname))
				stop(msg, call. = FALSE)
			}
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
		# Check the data has row and column names as appropriate
		verifyObservedNames(mxDataObject@observed, mxDataObject@means, mxDataObject@type, flatModel, modelname, "LISREL")
		#
		# Change *Matrix from the string name of the matrix to the object
		lxMatrix <- flatModel[[lxMatrix]]
		lyMatrix <- flatModel[[lyMatrix]]
		beMatrix <- flatModel[[beMatrix]]
		gaMatrix <- flatModel[[gaMatrix]]
		phMatrix <- flatModel[[phMatrix]]
		psMatrix <- flatModel[[psMatrix]]
		tdMatrix <- flatModel[[tdMatrix]]
		teMatrix <- flatModel[[teMatrix]]
		thMatrix <- flatModel[[thMatrix]]
		txMatrix <- flatModel[[txMatrix]]
		tyMatrix <- flatModel[[tyMatrix]]
		kaMatrix <- flatModel[[kaMatrix]]
		alMatrix <- flatModel[[alMatrix]]
		#
		# if LY is not null, then
		#	 check LY for dimnames and
		#	 check its means (TY, AL) and
		#	 check for TE, BE, PS
		if(!is.null(lyMatrix)){
			# Check LY for dimnames
			if (is.null(dimnames(lyMatrix))) {
				msg <- paste("The LY matrix of model",
					omxQuotes(modelname), "does not contain dimnames")
				stop(msg, call. = FALSE)
			}
			if (is.null(dimnames(lyMatrix)[[2]])) {
				msg <- paste("The LY matrix of model",
					omxQuotes(modelname), "does not contain colnames")
				stop(msg, call. = FALSE)
			}
			# Check its means (TY, AL) 
			if(!single.na(mxDataObject@means) || mxDataObject@type == "raw") {
				checkLISRELmeans(
					Lam=lyMatrix,
					ManMean=tyMatrix,
					LatMean=alMatrix,
					X=FALSE,
					modelname=modelname
				)
			}
			# Check for TE, BE, PS
			if(is.null(teMatrix)){
				msg <- paste("The TE matrix is absent but the LY",
					"matrix is present in model",
					omxQuotes(modelname))
				stop(msg, call. = FALSE)
			}
			if(is.null(beMatrix)){
				msg <- paste("The BE matrix is absent but the LY",
					"matrix is present in model",
					omxQuotes(modelname))
				stop(msg, call. = FALSE)
			}
			if(is.null(psMatrix)){
				msg <- paste("The PS matrix is absent but the LY",
					"matrix is present in model",
					omxQuotes(modelname))
				stop(msg, call. = FALSE)
			}
		}
		#
		# if LX is non-null, then
		#	 check LX for dimnames
		#	 check its means (TX, KA) and
		#	 check for TD, PH
		if(!is.null(lxMatrix)){
			# Check LX for dimnames
			if (is.null(dimnames(lxMatrix))) {
				msg <- paste("The LX matrix of model",
					omxQuotes(modelname), "does not contain dimnames")
				stop(msg, call. = FALSE)
			}
			if (is.null(dimnames(lxMatrix)[[2]])) {
				msg <- paste("The LX matrix of model",
					omxQuotes(modelname), "does not contain colnames")
				stop(msg, call. = FALSE)
			}
			# Check its means (TX, KA)
			if(!single.na(mxDataObject@means) || mxDataObject@type == "raw") {
				checkLISRELmeans(
					Lam=lxMatrix,
					ManMean=txMatrix,
					LatMean=kaMatrix,
					X=TRUE,
					modelname=modelname
				)
			}
			# Check for TD, PH
			if(is.null(tdMatrix)){
				msg <- paste("The TD matrix is absent but the LX",
					"matrix is present in model",
					omxQuotes(modelname))
				stop(msg, call. = FALSE)
			}
			if(is.null(phMatrix)){
				msg <- paste("The PH matrix is absent but the LX",
					"matrix is present in model",
					omxQuotes(modelname))
				stop(msg, call. = FALSE)
			}
		}
		#
		# if both LX and LY are not null
		#	 must have TH, GA
		if(!is.null(lxMatrix) && !is.null(lyMatrix)){
			if(is.null(thMatrix)){
				msg <- paste("The TH matrix is absent but the LY and LX",
					"matrices are present in model",
					omxQuotes(modelname))
				stop(msg, call. = FALSE)
			}
			if(is.null(gaMatrix)){
				msg <- paste("The GA matrix is absent but the LY and LX",
					"matrices are present in model",
					omxQuotes(modelname))
				stop(msg, call. = FALSE)
			}
		}
		if(is.null(lxMatrix) && (!is.null(phMatrix) || !is.null(tdMatrix))) {
			msg <- paste("Some, but not all of measurement and structural matrices are missing.",
				"The LX matrix is absent but the PH or TD",
				"matrices are present in model",
				omxQuotes(modelname))
			stop(msg, call. = FALSE)
		}
		if(is.null(lyMatrix) && (!is.null(psMatrix) || !is.null(teMatrix) || !is.null(beMatrix))) {
			msg <- paste("Some, but not all of measurement and structural matrices are missing.",
				"The LY matrix is absent but at least one of the PS, TE, and BE",
				"matrices are present in model",
				omxQuotes(modelname))
			stop(msg, call. = FALSE)
		}
		# Check if both LX and LY are missing
		if(is.null(lyMatrix) && is.null(lxMatrix)){
			msg <- paste("The model", omxQuotes(modelname), "must have at least",
				"one of the 'LX' matrix or the 'LY' matrix, but both appear to be missing.")
			stop(msg, call. = FALSE)
		}
		#
		# Raw data error checking
		#  Set the canonical order of observed variable names.
		translatedNames <- c(dimnames(lyMatrix)[[1]], dimnames(lxMatrix)[[1]]) #fMatrixTranslateNames(fMatrix, modelname) #Rearrange the rownames of F to match the order of the columns
		.Object@depth <- generateLISRELDepth(flatModel, beMatrix2, model@options) #Find out how many iterations of I + BE + BE^2 + ... are need until nilpotency.
		if (mxDataObject@type == 'raw') {
			threshName <- .Object@thresholds
			checkNumberOrdinalColumns(mxDataObject)
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
		} else {# Non-Raw data checking
			.Object@thresholds <- as.integer(NA)
			# Check the observed covariance matrix is separated into endo and exo blocks
			if (!identical(translatedNames, rownames(mxDataObject@observed))) {
				msg <- paste("The names of the manifest",
					"variables in the LY and LX matrices of model",
					omxQuotes(modelname), "do not match the",
					"dimnames of the observed covariance matrix",
					"or they are in the wrong order.")
				stop(msg, call. = FALSE)
			}
		}
		return(.Object)
	}
)


#--------------------------------------------------------------------
# **DONE**
setMethod("genericExpDependencies", signature("MxExpectationLISREL"),
	function(.Object, dependencies) {
	sources <- c(.Object@LX, .Object@LY, .Object@BE, .Object@GA, 
		.Object@PH, .Object@PS, .Object@TD, .Object@TE, 
		.Object@TH, .Object@TX, .Object@TY, .Object@KA, 
		.Object@AL, .Object@thresholds)
	sources <- sources[!is.na(sources)]
	dependencies <- imxAddDependency(sources, .Object@name, dependencies)
	return(dependencies)
	}
)


#--------------------------------------------------------------------
# **DONE**
setMethod("genericExpRename", signature("MxExpectationLISREL"),
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


#--------------------------------------------------------------------
checkLISRELargument <- function(x, xname) {
	if (!(single.na(x) || typeof(x) == "character")) {
		msg <- paste("argument ", xname, " is not a string ",
			"(the name of the '", xname, "' matrix)", sep="")
		stop(msg)
	}
	if (is.na(x)) x <- as.integer(NA)
	return(x)
}


#--------------------------------------------------------------------
# **DONE**
mxExpectationLISREL <- function(LX=NA, LY=NA, BE=NA, GA=NA, PH=NA, PS=NA, TD=NA, TE=NA, TH=NA, TX = NA, TY = NA, KA = NA, AL = NA, dimnames = NA, thresholds = NA, threshnames = dimnames) {
	LX <- checkLISRELargument(LX, "LX")
	LY <- checkLISRELargument(LY, "LY")
	BE <- checkLISRELargument(BE, "BE")
	GA <- checkLISRELargument(GA, "GA")
	PH <- checkLISRELargument(PH, "PH")
	PS <- checkLISRELargument(PS, "PS")
	TD <- checkLISRELargument(TD, "TD")
	TE <- checkLISRELargument(TE, "TE")
	TH <- checkLISRELargument(TH, "TH")
	TX <- checkLISRELargument(TX, "TX")
	TY <- checkLISRELargument(TY, "TY")
	KA <- checkLISRELargument(KA, "KA")
	AL <- checkLISRELargument(AL, "AL")
	
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
	return(new("MxExpectationLISREL", LX, LY, BE, GA, PH, PS, TD, TE, TH, TX, TY, KA, AL, dimnames, thresholds, threshnames))
}


#--------------------------------------------------------------------
# **DONE**
displayExpectationLISREL <- function(expectation) {
	cat("MxExpectationLISREL", omxQuotes(expectation@name), '\n')
	cat("$LX :", omxQuotes(expectation@LX), '\n')
	cat("$LY :", omxQuotes(expectation@LY), '\n')
	cat("$BE :", omxQuotes(expectation@BE), '\n')
	cat("$GA :", omxQuotes(expectation@GA), '\n')
	cat("$PH :", omxQuotes(expectation@PH), '\n')
	cat("$PS :", omxQuotes(expectation@PS), '\n')
	cat("$TD :", omxQuotes(expectation@TD), '\n')
	cat("$TE :", omxQuotes(expectation@TE), '\n')
	cat("$TH :", omxQuotes(expectation@TH), '\n')
	if (is.na(expectation@TX)) {
		cat("$TX :", expectation@TX, '\n')
	} else {
		cat("$TX :", omxQuotes(expectation@TX), '\n')
	}
	if (is.na(expectation@TY)) {
		cat("$TY :", expectation@TY, '\n')
	} else {
		cat("$TY :", omxQuotes(expectation@TY), '\n')
	}
	if (is.na(expectation@KA)) {
		cat("$KA :", expectation@KA, '\n')
	} else {
		cat("$KA :", omxQuotes(expectation@KA), '\n')
	}
	if (is.na(expectation@AL)) {
		cat("$AL :", expectation@AL, '\n')
	} else {
		cat("$AL :", omxQuotes(expectation@AL), '\n')
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


#--------------------------------------------------------------------
# **DONE**
setMethod("print", "MxExpectationLISREL", function(x,...) { 
	displayExpectationLISREL(x) 
})


#--------------------------------------------------------------------
# **DONE**
setMethod("show", "MxExpectationLISREL", function(object) { 
	displayExpectationLISREL(object) 
})


#------------------------------------------------------------------------------
setMethod("genericGetExpected", signature("MxExpectationLISREL"),
	  function(.Object, model, what, defvar.row=1) {
		  ret <- list()
		  LXname <- .Object@LX
		  LYname <- .Object@LY
		  BEname <- .Object@BE
		  GAname <- .Object@GA
		  PHname <- .Object@PH
		  PSname <- .Object@PS
		  TDname <- .Object@TD
		  TEname <- .Object@TE
		  THname <- .Object@TH
		  TXname <- .Object@TX
		  TYname <- .Object@TY
		  KAname <- .Object@KA
		  ALname <- .Object@AL
		  hasX <- !single.na(LXname)
		  hasY <- !single.na(LYname)
		  if(hasX){
			  LX <- mxEvalByName(LXname, model, compute=TRUE, defvar.row=defvar.row)
			  PH <- mxEvalByName(PHname, model, compute=TRUE, defvar.row=defvar.row)
			  TD <- mxEvalByName(TDname, model, compute=TRUE, defvar.row=defvar.row)
		  } else {
			  LX <- matrix( , 0, 0)
			  PH <- matrix( , 0, 0)
			  TD <- matrix( , 0, 0)
		  }
		  if(hasY){
			  LY <- mxEvalByName(LYname, model, compute=TRUE, defvar.row=defvar.row)
			  BE <- mxEvalByName(BEname, model, compute=TRUE, defvar.row=defvar.row)
			  PS <- mxEvalByName(PSname, model, compute=TRUE, defvar.row=defvar.row)
			  TE <- mxEvalByName(TEname, model, compute=TRUE, defvar.row=defvar.row)
			  AL <- mxEvalByName(ALname, model, compute=TRUE, defvar.row=defvar.row)
			  I <- diag(1, nrow=nrow(BE))
			  A <- LY %*% solve(I-BE)
		  } else {
			  LY <- matrix( , 0, 0)
			  BE <- matrix( , 0, 0)
			  PS <- matrix( , 0, 0)
			  TE <- matrix( , 0, 0)
			  AL <- matrix( , 0, 1)
			  A <- matrix( , 0, 0)
		  }
		  if(hasX & hasY){
			  GA <- mxEvalByName(GAname, model, compute=TRUE, defvar.row=defvar.row)
			  TH <- mxEvalByName(THname, model, compute=TRUE, defvar.row=defvar.row)
		  } else {
			  GA <- matrix( , nrow=ncol(LY), ncol=ncol(LX))
			  TH <- matrix( , nrow=nrow(LX), ncol=nrow(LY))
		  }
		  if ('covariance' %in% what) {
			  endoBlock <- A %*% (GA %*% PH %*% t(GA) + PS) %*% t(A) + TE
			  exoBlock <- LX %*% PH %*% t(LX) + TD
			  exenBlock <- LX %*% PH %*% t(GA) %*% t(A) + TH
			  cov <- rbind(cbind(endoBlock, t(exenBlock)),
				       cbind(exenBlock, exoBlock))
			  ret[['covariance']] <- cov
		  }
		  if ('means' %in% what) {
			  if(single.na(TXname) & single.na(TYname)){
				  mean <- matrix( , 0, 0)
			  } else {
					if(hasX & single.na(TXname)){
						stop("Model has exogenous variables but not exogenous means.")
					}
					if(hasY & single.na(TYname)){
						stop("Model has endogenous variables but not endogenous means.")
					}
					endoMean <- NULL
					exoMean <- NULL
					if(hasY){
						KA <- mxEvalByName(KAname, model, compute=TRUE, defvar.row=defvar.row)
						TY <- mxEvalByName(TYname, model, compute=TRUE, defvar.row=defvar.row)
						endoMean <- TY + A %*% (AL + GA %*% KA)
					}
					if(hasX){
						TX <- mxEvalByName(TXname, model, compute=TRUE, defvar.row=defvar.row)
						KA <- mxEvalByName(KAname, model, compute=TRUE, defvar.row=defvar.row)
						exoMean <- TX + LX %*% KA
					}
					mean <- rbind(endoMean, exoMean)
			  }
			  ret[['means']] <- mean
		  }
		  if ('thresholds' %in% what) {
			  thrname <- .Object@thresholds
			  if(!single.na(thrname)){
				  thr <- mxEvalByName(thrname, model, compute=TRUE, defvar.row=defvar.row)
			  } else {thr <- matrix( , 0, 0)}
			  ret[['thresholds']] <- thr
		  }
		  ret
})


setMethod("genericExpGetPrecision", "MxExpectationLISREL",
	function(.Object) {
		if(!single.na(.Object@thresholds)) {
			return(list(stepSize=0.1, iterations=3L))
		} else {
			callNextMethod();
		}
})



#------------------------------------------------------------------------------
setMethod("genericGenerateData", signature("MxExpectationLISREL"),
	function(.Object, model, nrows) {
		return(generateNormalData(model, nrows))
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
generateLISRELDepth <- function(flatModel, aMatrixName, modeloptions) {
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
