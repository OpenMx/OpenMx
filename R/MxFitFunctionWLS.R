#
#   Copyright 2007-2013 The OpenMx Project
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

# TODO: Add conformance and dimname checking for weight matrix

#------------------------------------------------------------------------------
# Revision History
#  Fri 08 Mar 2013 15:35:29 Central Standard Time -- Michael Hunter copied file from 
#    OpenMx\branches\dependency-tracking\R\WLSObjective.R
#    I think Tim Brick wrote the original.
#    Michael Hunter edited the file to use fitfunctions
#  


#------------------------------------------------------------------------------

# **DONE**
setClass(Class = "MxFitFunctionWLS",
	representation = representation(weights = "MxCharOrNumber"),
	contains = "MxBaseFitFunction")

# **DONE**
setMethod("initialize", "MxFitFunctionWLS",
	function(.Object, weights, name = 'fitfunction') {
		.Object@name <- name
		.Object@weights <- weights
		return(.Object)
	}
)

# **DONE**
setMethod("qualifyNames", signature("MxFitFunctionWLS"), 
	function(.Object, modelname, namespace) {
		.Object@name <- imxIdentifier(modelname, .Object@name)
		.Object@weights <- imxConvertIdentifier(.Object@weights,
			modelname, namespace)
		return(.Object)
})

# **DONE**
setMethod("genericFitConvertEntities", "MxFitFunctionWLS", 
	function(.Object, flatModel, namespace, labelsData) {
		name <- .Object@name
		modelname <- imxReverseIdentifier(model, .Object@name)[[1]]
		expectName <- paste(modelname, "expectation", sep=".")
		
		expectation <- flatModel@expectations[[expectName]]
		dataname <- expectation@data
		
		if (flatModel@datasets[[dataname]]@type != 'raw') {
			if (.Object@vector) {
				modelname <- getModelName(.Object)
				msg <- paste("The WLS fit function",
					"in model", omxQuotes(modelname),
					"does not have raw data")
				stop(msg, call.=FALSE)
			}
		}
		return(flatModel)
})

# **DONE**
setMethod("genericFitFunConvert", "MxFitFunctionWLS", 
	function(.Object, flatModel, model, labelsData, defVars, dependencies) {
		name <- .Object@name
		modelname <- imxReverseIdentifier(model, .Object@name)[[1]]
		expectName <- paste(modelname, "expectation", sep=".")
		if (expectName %in% names(flatModel@expectations)) {
			expectIndex <- imxLocateIndex(flatModel, expectName, name)
		} else {
			expectIndex <- as.integer(NA)
		}
		.Object@expectation <- expectIndex
		return(.Object)
})

# **DONE**
setMethod("genericFitInitialMatrix", "MxFitFunctionWLS",
	function(.Object, flatModel) {
		flatFitFunction <- flatModel@fitfunctions[[.Object@name]]
		modelname <- imxReverseIdentifier(flatModel, flatFitFunction@name)[[1]]
		expectationName <- paste(modelname, "expectation", sep = ".")
		expectation <- flatModel@expectations[[expectationName]]
		if (is.null(expectation)) {
			msg <- paste("The WLS fit function",
			"has a missing expectation in the model",
			omxQuotes(modelname))
			stop(msg, call.=FALSE)
		}
		if (is.na(expectation@data)) {
			msg <- paste("The WLS fit function has",
			"an expectation function with no data in the model",
			omxQuotes(modelname))
			stop(msg, call.=FALSE)
		}
		mxDataObject <- flatModel@datasets[[expectation@data]]
		if (mxDataObject@type != 'raw') {
			msg <- paste("The dataset associated with the WLS expectation function", 
				"in model", omxQuotes(modelname), "is not raw data.")
			stop(msg, call.=FALSE)
		}
		cols <- ncol(mxDataObject@observed)
		return(matrix(as.double(NA), cols, cols))
})

# **DONE**
mxFitFunctionWLS <- function(weights="ULS") {
	if (is.na(weights)) weights <- as.character(NA)
	if (!(weights %in% c("ULS", "DWLS", "WLS"))){
		stop("weights must be either 'ULS', 'DLS' or 'WLS'")
		}
	return(new("MxFitFunctionWLS", weights))
}

# **DONE**
displayMxFitFunctionWLS <- function(fitfunction) {
	cat("MxFitFunctionWLS", omxQuotes(fitfunction@name), '\n')
	if (single.na(fitfunction@weights)) {
		cat("@weights : NA \n")
	} else {
		cat("@weights :", omxQuotes(fitfunction@weights), '\n')
	}
	if (length(fitfunction@result) == 0) {
		cat("@result: (not yet computed) ")
	} else {
		cat("@result:\n")
	}
	print(fitfunction@result)
	invisible(fitfunction)
}

# **DONE**
setMethod("print", "MxFitFunctionWLS", function(x, ...) { 
	displayMxFitFunctionWLS(x) 
})

setMethod("show", "MxFitFunctionWLS", function(object) { 
	displayMxFitFunctionWLS(object) 
})

