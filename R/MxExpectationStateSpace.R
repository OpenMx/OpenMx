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
# Filename: MxExpectationStateSpace.R
# Date: 2012.11.14
# Purpose: Define classes and methods for the state space model (SSM)
#  expectations.
#--------------------------------------------------------------------

#--------------------------------------------------------------------
# Revision History
#   Wed Nov 14 13:34:01 Central Standard Time 2012 -- Michael Hunter created file
#   Sat Nov 17 16:14:56 Central Standard Time 2012 -- Michael Hunter change names to ExpectationStateSpace instead of ExpectationSSM
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
	function(.Object, A, B, C, D, Q, R, dims, thresholds, threshnames,
		data = as.integer(NA), name = 'expectation') {
		.Object@name <- name
		.Object@A <- A
		.Object@B <- B
		.Object@C <- C
		.Object@D <- D
		.Object@Q <- Q
		.Object@R <- R
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
		if(is.na(.Object@data)) {
			msg <- paste("The SSM expectation function",
				"does not have a dataset associated with it in model",
				omxQuotes(model@name))
			stop(msg, call.=FALSE)
		}
		flatModel <- updateThresholdDimnames(.Object, flatModel, labelsData)
		
		return(flatModel)
	}
)


#--------------------------------------------------------------------
# **DONE**
setMethod("genericExpFunNamespace", signature("MxExpectationStateSpace"), 
	function(.Object, modelname, namespace) {
		.Object@name <- imxIdentifier(modelname, .Object@name)
		.Object@A <- imxConvertIdentifier(.Object@A, modelname, namespace)
		.Object@B <- imxConvertIdentifier(.Object@B, modelname, namespace)
		.Object@C <- imxConvertIdentifier(.Object@C, modelname, namespace)
		.Object@D <- imxConvertIdentifier(.Object@D, modelname, namespace)
		.Object@Q <- imxConvertIdentifier(.Object@Q, modelname, namespace)
		.Object@R <- imxConvertIdentifier(.Object@R, modelname, namespace)
		.Object@data <- imxConvertIdentifier(.Object@data, modelname, namespace)
		.Object@thresholds <- sapply(.Object@thresholds, imxConvertIdentifier, modelname, namespace)
		return(.Object)
	}
)


#--------------------------------------------------------------------
# TODO: Add lots of error checking.
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
		data <- .Object@data
		if(is.na(data)) {
			msg <- paste("The SSM expectation function",
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
		.Object@data <- as.integer(imxLocateIndex(flatModel, data, name))
		return(.Object)
	}
)


#--------------------------------------------------------------------
# **DONE**
setMethod("genericExpDependencies", signature("MxExpectationStateSpace"),
	function(.Object, dependencies) {
		sources <- c(.Object@A, .Object@B, .Object@C, .Object@D, .Object@Q, .Object@R, .Object@thresholds)
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
		.Object@data <- renameReference(.Object@data, oldname, newname)
		.Object@thresholds <- sapply(.Object@thresholds, renameReference, oldname, newname)		
		return(.Object)
	}
)


#--------------------------------------------------------------------
# Note: this turns off data sorting for the State Space expectation
setMethod("genericExpAddEntities", "MxExpectationStateSpace",
        function(.Object, job, flatJob, labelsData) {
                key <- "No Sort Data"
                value <- getModelName(.Object)
                job <- mxOption(job, key, value)
                return(job)
        }
)


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
imxExpectationStateSpace <- function(A, B, C, D, Q, R, dimnames = NA, thresholds = NA, threshnames = dimnames){
	A <- checkLISRELargument(A, "A")
	B <- checkLISRELargument(B, "B")
	C <- checkLISRELargument(C, "C")
	D <- checkLISRELargument(D, "D")
	Q <- checkLISRELargument(Q, "Q")
	R <- checkLISRELargument(R, "R")
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
	return(new("MxExpectationStateSpace", A, B, C, D, Q, R, dimnames, thresholds, threshnames))
}


#--------------------------------------------------------------------
# TODO: Add expected mean and cov printouts
displayMxExpectationStateSpace <- function(expectation) {
	cat("MxExpectationStateSpace", omxQuotes(expectation@name), '\n')
	cat("@A :", omxQuotes(expectation@A), '\n')
	cat("@B :", omxQuotes(expectation@B), '\n')
	cat("@C :", omxQuotes(expectation@C), '\n')
	cat("@D :", omxQuotes(expectation@D), '\n')
	cat("@Q :", omxQuotes(expectation@Q), '\n')
	cat("@R :", omxQuotes(expectation@R), '\n')
	if (single.na(expectation@dims)) {
		cat("@dims : NA \n")
	} else {
		cat("@dims :", omxQuotes(expectation@dims), '\n')
	}		
	if (single.na(expectation@thresholds)) {
		cat("@thresholds : NA \n")
	} else {
		cat("@thresholds :", omxQuotes(expectation@thresholds), '\n')
	}
	cat("@info$likelihoods: ", length(expectation@info$likelihoods) > 0, '\n')
	invisible(expectation)
}

setMethod("print", "MxExpectationStateSpace", function(x,...) { 
	displayMxExpectationStateSpace(x) 
})

setMethod("show", "MxExpectationStateSpace", function(object) { 
	displayMxExpectationStateSpace(object) 
})


