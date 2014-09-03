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


#--------------------------------------------------------------------
# Author: Michael D. Hunter
# Filename: MxExpectationStateSpace.R
# Date: 2012.11.14
# Purpose: Define classes and methods for the state space model (SSM)
#  expectations.
#--------------------------------------------------------------------

#--------------------------------------------------------------------
# Revision History
#   Wed Nov 14 13:34:01 Central Standard Time 2012 -- Michael Hunter created file
#   Sat Nov 17 16:14:56 Central Standard Time 2012 -- Michael Hunter change names to ExpectationStateSpace instead of ExpectationSSM
#   Wed Sep 18 23:09:34 Central Daylight Time 2013 -- Michael Hunter changed name of arguments to x0 and P0
#   

#--------------------------------------------------------------------
# **DONE**
setClass(Class = "MxExpectationStateSpace",
	representation = representation(
		A = "MxCharOrNumber",
		B = "MxCharOrNumber",
		C = "MxCharOrNumber",
		D = "MxCharOrNumber",
		Q = "MxCharOrNumber",
		R = "MxCharOrNumber",
		x0 = "MxCharOrNumber",
		P0 = "MxCharOrNumber",
		u = "MxCharOrNumber",
		thresholds = "MxCharOrNumber",
		dims = "character",
		definitionVars = "list",
		dataColumns = "numeric",
		thresholdColumns = "numeric",
		thresholdLevels = "numeric",
		threshnames = "character"),
	contains = "MxBaseExpectation")


#--------------------------------------------------------------------
# **DONE**
setMethod("initialize", "MxExpectationStateSpace",
	function(.Object, A, B, C, D, Q, R, x0, P0, u, dims, thresholds, threshnames,
		data = as.integer(NA), name = 'expectation') {
		.Object@name <- name
		.Object@A <- A
		.Object@B <- B
		.Object@C <- C
		.Object@D <- D
		.Object@Q <- Q
		.Object@R <- R
		.Object@x0 <- x0
		.Object@P0 <- P0
		.Object@u <- u
		.Object@data <- data
		.Object@dims <- dims
		.Object@thresholds <- thresholds
		.Object@definitionVars <- list()
		.Object@threshnames <- threshnames
		return(.Object)
	}
)


#--------------------------------------------------------------------
# TODO: Ask Spiegel and Brick what this is really supposed to do.
setMethod("genericExpConvertEntities", "MxExpectationStateSpace",
	function(.Object, flatModel, namespace, labelsData) {
		cache <- new.env(parent = emptyenv())
		modelname <- getModelName(.Object)
		if(is.na(.Object@data)) {
			msg <- paste("The SSM expectation function",
				"does not have a dataset associated with it in model",
				omxQuotes(modelname))
			stop(msg, call.=FALSE)
		}
		
		tuple <- evaluateMxObject(.Object@A, flatModel, labelsData, cache)
		Amatrix <- tuple[[1]]
		cache <- tuple[[2]]
		tuple <- evaluateMxObject(.Object@B, flatModel, labelsData, cache)
		Bmatrix <- tuple[[1]]
		cache <- tuple[[2]]
		tuple <- evaluateMxObject(.Object@C, flatModel, labelsData, cache)
		Cmatrix <- tuple[[1]]
		cache <- tuple[[2]]
		tuple <- evaluateMxObject(.Object@D, flatModel, labelsData, cache)
		Dmatrix <- tuple[[1]]
		cache <- tuple[[2]]
		tuple <- evaluateMxObject(.Object@Q, flatModel, labelsData, cache)
		Qmatrix <- tuple[[1]]
		cache <- tuple[[2]]
		tuple <- evaluateMxObject(.Object@R, flatModel, labelsData, cache)
		Rmatrix <- tuple[[1]]
		cache <- tuple[[2]]
		tuple <- evaluateMxObject(.Object@x0, flatModel, labelsData, cache)
		Xmatrix <- tuple[[1]]
		cache <- tuple[[2]]
		tuple <- evaluateMxObject(.Object@P0, flatModel, labelsData, cache)
		Pmatrix <- tuple[[1]]
		cache <- tuple[[2]]
		tuple <- evaluateMxObject(.Object@u, flatModel, labelsData, cache)
		Umatrix <- tuple[[1]]
		cache <- tuple[[2]]
		# Conformability checking
		ldim <- ncol(Cmatrix)
		mdim <- nrow(Cmatrix)
		checkSSMConformable(Amatrix, ldim, ldim, 'A', omxQuotes(modelname))
		checkSSMConformable(Cmatrix, mdim, ldim, 'C', omxQuotes(modelname))
		checkSSMConformable(Qmatrix, ldim, ldim, 'Q', omxQuotes(modelname))
		checkSSMConformable(Rmatrix, mdim, mdim, 'R', omxQuotes(modelname))
		checkSSMConformable(Xmatrix, ldim, 1, 'x0', omxQuotes(modelname))
		checkSSMConformable(Pmatrix, ldim, ldim, 'P0', omxQuotes(modelname))
		#
		# If any of B, D, u are not missing, then
		#  1.  All of B, D, u must not be missing
		#  2.  colnames of B and D must match rownames of u?
		# 3.   Conformability checks
		#  ldim <- ncol(Cmatrix)
		#  mdim <- nrow(Cmatrix)
		#  udim <- ncol(Bmatrix)
		#  A: ldim x ldim
		#  B: ldim x udim
		#  C: mdim x ldim
		#  D: mdim x udim
		#  Q: ldim x ldim
		#  R: mdim x mdim
		#  x: ldim x 1
		#  P: ldim x ldim
		#  u: udim x 1
		if(!is.null(Bmatrix) || !is.null(Dmatrix) || !is.null(Umatrix)){
			udim <- ncol(Bmatrix)
			checkSSMConformable(Bmatrix, ldim, udim, 'B', omxQuotes(modelname))
			checkSSMConformable(Dmatrix, mdim, udim, 'D', omxQuotes(modelname))
			checkSSMConformable(Umatrix, udim, 1, 'u', omxQuotes(modelname))
		}
		flatModel <- updateSSMdimnames(.Object, flatModel)
		flatModel <- updateThresholdDimnames(.Object, flatModel, labelsData)
		
		return(flatModel)
	}
)

updateSSMdimnames <- function(flatExpectation, flatJob){
	cMatrixName <- flatExpectation@C
	cMatrix <- flatJob[[cMatrixName]]
	if (is.null(cMatrix)) {
		modelname <- getModelName(flatExpectation)
		stop(paste("Unknown C matrix name", 
			omxQuotes(simplifyName(cMatrixName, modelname)),
			"detected in the State Space expectation function",
			"of model", omxQuotes(modelname)), call. = FALSE)
	}
	dims <- flatExpectation@dims
	if (!is.null(dimnames(cMatrix)) && !single.na(dims) && 
		!identical(dimnames(cMatrix)[[1]], dims)) {
		modelname <- getModelName(flatExpectation)
		msg <- paste("The C matrix associated",
			"with the State Space expectation function in model", 
			omxQuotes(modelname), "contains dimnames and",
			"the expectation function has specified different dimnames")
		stop(msg, call.=FALSE)		
	}
	if (is.null(dimnames(cMatrix)) && !single.na(dims)) {
		dimnames(flatJob[[cMatrixName]]) <- list(dims, c())
	}
	
	return(flatJob)
}


#--------------------------------------------------------------------
# **DONE**
setMethod("qualifyNames", signature("MxExpectationStateSpace"), 
	function(.Object, modelname, namespace) {
		.Object@name <- imxIdentifier(modelname, .Object@name)
		.Object@A <- imxConvertIdentifier(.Object@A, modelname, namespace)
		.Object@B <- imxConvertIdentifier(.Object@B, modelname, namespace)
		.Object@C <- imxConvertIdentifier(.Object@C, modelname, namespace)
		.Object@D <- imxConvertIdentifier(.Object@D, modelname, namespace)
		.Object@Q <- imxConvertIdentifier(.Object@Q, modelname, namespace)
		.Object@R <- imxConvertIdentifier(.Object@R, modelname, namespace)
		.Object@x0 <- imxConvertIdentifier(.Object@x0, modelname, namespace)
		.Object@P0 <- imxConvertIdentifier(.Object@P0, modelname, namespace)
		.Object@u <- imxConvertIdentifier(.Object@u, modelname, namespace)
		.Object@data <- imxConvertIdentifier(.Object@data, modelname, namespace)
		.Object@thresholds <- sapply(.Object@thresholds, imxConvertIdentifier, modelname, namespace)
		return(.Object)
	}
)




#--------------------------------------------------------------------

checkSSMNotMissing <- function(matrixobj, matrixname, modelname){
	if(is.null(matrixobj)){
		msg <- paste("The", matrixname, "matrix of the state space expectation in model",
			modelname, "is missing.")
		stop(msg, call. = FALSE)
	}
}

checkSSMConformable <- function(mat, rows, cols, matname, modname){
	if( nrow(mat) != rows || ncol(mat) != cols ){
		msg <- paste("The ", matname, " matrix is not the correct size",
			" in the state space expectation of model ", modname,
			".  It is ", nrow(mat), " by ", ncol(mat), " and should be ",
			rows, " by ", cols, ".", sep="")
		stop(msg, call. = FALSE)
	}
}

# TODO: Allow subsets of the matrices to be specified
#  by filling in default matrices.
setMethod("genericExpFunConvert", signature("MxExpectationStateSpace"), 
	function(.Object, flatModel, model, labelsData, defVars, dependencies) {
		modelname <- imxReverseIdentifier(model, .Object@name)[[1]]	
		name <- .Object@name
		aMatrix <- .Object@A
		bMatrix <- .Object@B
		cMatrix <- .Object@C
		dMatrix <- .Object@D
		qMatrix <- .Object@Q
		rMatrix <- .Object@R
		xMatrix <- .Object@x0
		pMatrix <- .Object@P0
		uMatrix <- .Object@u
		data <- .Object@data
		if(is.na(data)) {
			msg <- paste("The state space expectation function",
				"does not have a dataset associated with it in model",
				omxQuotes(modelname))
			stop(msg, call. = FALSE)
		}
		mxDataObject <- flatModel@datasets[[.Object@data]]
		checkNumericData(mxDataObject)
		.Object@A <- imxLocateIndex(flatModel, aMatrix, name)
		.Object@B <- imxLocateIndex(flatModel, bMatrix, name)
		.Object@C <- imxLocateIndex(flatModel, cMatrix, name)
		.Object@D <- imxLocateIndex(flatModel, dMatrix, name)
		.Object@Q <- imxLocateIndex(flatModel, qMatrix, name)
		.Object@R <- imxLocateIndex(flatModel, rMatrix, name)
		.Object@x0 <- imxLocateIndex(flatModel, xMatrix, name)
		.Object@P0 <- imxLocateIndex(flatModel, pMatrix, name)
		.Object@u <- imxLocateIndex(flatModel, uMatrix, name)
		.Object@data <- as.integer(imxLocateIndex(flatModel, data, name))
		#
		# Check the data has row and column names as appropriate
		verifyObservedNames(mxDataObject@observed, mxDataObject@means, mxDataObject@type, flatModel, modelname, "State Space")
		#
		# Change *Matrix from the string name of the matrix to the object
		aMatrix <- flatModel[[aMatrix]]
		bMatrix <- flatModel[[bMatrix]]
		cMatrix <- flatModel[[cMatrix]]
		dMatrix <- flatModel[[dMatrix]]
		qMatrix <- flatModel[[qMatrix]]
		rMatrix <- flatModel[[rMatrix]]
		xMatrix <- flatModel[[xMatrix]]
		pMatrix <- flatModel[[pMatrix]]
		uMatrix <- flatModel[[uMatrix]]
		#
		# Check for missing matrices or combinations of matrices
		checkSSMNotMissing(cMatrix, 'C', omxQuotes(modelname))
		# Check C for dimnames
		if (is.null(dimnames(cMatrix))) {
			msg <- paste("The C matrix of model",
				omxQuotes(modelname), "does not contain dimnames")
			stop(msg, call. = FALSE)
		}
		if (is.null(dimnames(cMatrix)[[2]])) {
			msg <- paste("The C matrix of model",
				omxQuotes(modelname), "does not contain colnames")
			stop(msg, call. = FALSE)
		}
		checkSSMNotMissing(aMatrix, 'A', omxQuotes(modelname))
		checkSSMNotMissing(qMatrix, 'Q', omxQuotes(modelname))
		checkSSMNotMissing(rMatrix, 'R', omxQuotes(modelname))
		checkSSMNotMissing(xMatrix, 'x0', omxQuotes(modelname))
		checkSSMNotMissing(pMatrix, 'P0', omxQuotes(modelname))
		if(!is.null(bMatrix) || !is.null(dMatrix) || !is.null(uMatrix)){
			checkSSMNotMissing(bMatrix, 'B', omxQuotes(modelname))
			checkSSMNotMissing(dMatrix, 'D', omxQuotes(modelname))
			checkSSMNotMissing(uMatrix, 'u', omxQuotes(modelname))
		}
		translatedNames <- c(dimnames(cMatrix)[[1]])
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
		}
		return(.Object)
	}
)


#--------------------------------------------------------------------
# **DONE**
setMethod("genericExpDependencies", signature("MxExpectationStateSpace"),
	function(.Object, dependencies) {
		sources <- c(.Object@A, .Object@B, .Object@C, .Object@D, .Object@Q, .Object@R, .Object@x0, .Object@P0, .Object@u, .Object@thresholds)
		sources <- sources[!is.na(sources)]
		dependencies <- imxAddDependency(sources, .Object@name, dependencies)
		return(dependencies)
	}
)


#--------------------------------------------------------------------
# **DONE**
setMethod("genericExpRename", signature("MxExpectationStateSpace"),
	function(.Object, oldname, newname) {
		.Object@A <- renameReference(.Object@A, oldname, newname)
		.Object@B <- renameReference(.Object@B, oldname, newname)
		.Object@C <- renameReference(.Object@C, oldname, newname)
		.Object@D <- renameReference(.Object@D, oldname, newname)
		.Object@Q <- renameReference(.Object@Q, oldname, newname)
		.Object@R <- renameReference(.Object@R, oldname, newname)
		.Object@x0 <- renameReference(.Object@x0, oldname, newname)
		.Object@P0 <- renameReference(.Object@P0, oldname, newname)
		.Object@u <- renameReference(.Object@u, oldname, newname)
		.Object@data <- renameReference(.Object@data, oldname, newname)
		.Object@thresholds <- sapply(.Object@thresholds, renameReference, oldname, newname)		
		return(.Object)
	}
)


#--------------------------------------------------------------------
# Note: this turns off data sorting for the State Space expectation
setMethod("genericExpAddEntities", "MxExpectationStateSpace",
        function(.Object, job, flatJob, labelsData) {
                #TODO figure out how to handle situation where submodel with state space expectation
                # inherits its data from parent model.
                key <- "No Sort Data"
                #value <- c(job@options[[key]], getModelName(.Object)) #just add the model with the SSM exp to no sort
                value <- unique(c(job@options[[key]], getAllModelNames(job))) # add every model in the whole tree to to no sort
                # This is the nuclear option: whenever any model anywhere in the model tree has a SSMexp, don't sort any data.
                job <- mxOption(job, key, value)
                
                # Run state space models single threaded
                key <- "Number of Threads"
                job <- mxOption(job, key, 1)
                return(job)
        }
)

getAllModelNames <- function(model){
	ret <- getModelName(model)
	if(length(model@submodels) > 0){
		ret <- c(ret, unlist(lapply(model@submodels, getAllModelNames)))
	}
	return(ret)
}


#--------------------------------------------------------------------
checkSSMargument <- function(x, xname) {
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
mxExpectationStateSpace <- function(A, B, C, D, Q, R, x0, P0, u, dimnames = NA, thresholds = NA, threshnames = dimnames){
	A <- checkSSMargument(A, "A")
	B <- checkSSMargument(B, "B")
	C <- checkSSMargument(C, "C")
	D <- checkSSMargument(D, "D")
	Q <- checkSSMargument(Q, "Q")
	R <- checkSSMargument(R, "R")
	x0 <- checkSSMargument(x0, "x0")
	P0 <- checkSSMargument(P0, "P0")
	u <- checkSSMargument(u, "u")
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
	return(new("MxExpectationStateSpace", A, B, C, D, Q, R, x0, P0, u, dimnames, thresholds, threshnames))
}


#--------------------------------------------------------------------
# TODO: Add expected mean and cov printouts
displayMxExpectationStateSpace <- function(expectation) {
	cat("MxExpectationStateSpace", omxQuotes(expectation@name), '\n')
	cat("$A :", omxQuotes(expectation@A), '\n')
	cat("$B :", omxQuotes(expectation@B), '\n')
	cat("$C :", omxQuotes(expectation@C), '\n')
	cat("$D :", omxQuotes(expectation@D), '\n')
	cat("$Q :", omxQuotes(expectation@Q), '\n')
	cat("$R :", omxQuotes(expectation@R), '\n')
	cat("$x0 :", omxQuotes(expectation@x0), '\n')
	cat("$P0 :", omxQuotes(expectation@P0), '\n')
	cat("$u :", omxQuotes(expectation@u), '\n')
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

setMethod("print", "MxExpectationStateSpace", function(x,...) { 
	displayMxExpectationStateSpace(x) 
})

setMethod("show", "MxExpectationStateSpace", function(object) { 
	displayMxExpectationStateSpace(object) 
})


