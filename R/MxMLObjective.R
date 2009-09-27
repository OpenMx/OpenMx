#
#   Copyright 2007-2009 The OpenMx Project
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


setClass(Class = "MxMLObjective",
	representation = representation(
		covariance = "MxCharOrNumber",
		means = "MxCharOrNumber",
		thresholds = "character"),
	contains = "MxBaseObjective")

setMethod("initialize", "MxMLObjective",
	function(.Object, covariance, means, thresholds, 
		data = as.integer(NA), name = 'objective') {
		.Object@name <- name
		.Object@covariance <- covariance
		.Object@means <- means
		.Object@data <- data
		.Object@thresholds <- thresholds
		.Object@dependencies <- c('covariance', 'means', 'thresholds')
		return(.Object)
	}
)

setMethod("omxObjFunNamespace", signature("MxMLObjective"), 
	function(.Object, modelname, namespace) {
		.Object@name <- omxIdentifier(modelname, .Object@name)
		.Object@covariance <- omxConvertIdentifier(.Object@covariance, 
			modelname, namespace)
		.Object@means <- omxConvertIdentifier(.Object@means, 
			modelname, namespace)
		.Object@data <- omxConvertIdentifier(.Object@data, 
			modelname, namespace)
		.Object@thresholds <- omxConvertIdentifier(.Object@thresholds, 
			modelname, namespace)			
		return(.Object)
})

setMethod("omxObjFunConvert", signature("MxMLObjective"), 
	function(.Object, flatModel, model) {
		name <- .Object@name
		covariance <- .Object@covariance
		means <- .Object@means
		data <- .Object@data
		thresholds <- .Object@thresholds
		covarianceIndex <- omxLocateIndex(flatModel, covariance, name)
		if(is.na(data)) {
			msg <- paste("In model", omxQuotes(model@name),
				"the ML objective does not have a dataset specified")
			stop(msg, call. = FALSE)
		}
		mxDataObject <- flatModel@datasets[[.Object@data]]
		checkNumericData(mxDataObject)
		verifyExpectedNames(covariance, means, flatModel, "ML")
		meansIndex <- omxLocateIndex(flatModel, means, name)
		dIndex <- omxLocateIndex(flatModel, data, name)
		.Object@covariance <- covarianceIndex
		.Object@means <- meansIndex
		.Object@data <- dIndex
		return(.Object)
})


setMethod("omxObjModelConvert", "MxMLObjective",
	function(.Object, flatModel, model) {
		if(is.na(.Object@data)) {
			msg <- paste("The ML objective",
				"does not have a dataset associated with it in model",
				omxQuotes(model@name))
			stop(msg, call. = FALSE)
		}		
		if (flatModel@datasets[[.Object@data]]@type != 'raw') {
			return(model)
		}
		if (is.na(.Object@means)) {
			msg <- paste("In model", omxQuotes(model@name),
				"the ML objective has a raw dataset specified",
				"but no expected means vector")
			stop(msg, call. = FALSE)
		}
		objective <- mxFIMLObjective(.Object@covariance, .Object@means, .Object@thresholds)
		model@objective <- objective
		return(model)
	}
)


mxMLObjective <- function(covariance, means = NA, thresholds = NA) {
	if (missing(covariance) || typeof(covariance) != "character") {
		stop("Covariance argument is not a string (the name of the expected covariance matrix)")
	}
	if (!(single.na(means) || typeof(means) == "character")) {
		stop("Means argument is not a string (the name of the expected means matrix)")
	}
	if (is.na(means)) means <- as.integer(NA)
	if (is.na(thresholds)) thresholds <- as.character(NA)
	return(new("MxMLObjective", covariance, means, thresholds))
}

displayMLObjective <- function(objective) {
	cat("MxMLObjective", omxQuotes(objective@name), '\n')
	cat("@covariance :", omxQuotes(objective@covariance), '\n')
	cat("@means :", omxQuotes(objective@means), '\n')
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
	invisible(objective)
}


setMethod("print", "MxMLObjective", function(x,...) { 
	displayMLObjective(x) 
})

setMethod("show", "MxMLObjective", function(object) { 
	displayMLObjective(object) 
})
