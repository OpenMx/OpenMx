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
		definitionVars = "list"),
	contains = "MxBaseObjective")

setMethod("initialize", "MxMLObjective",
	function(.Object, covariance, means, data = as.numeric(NA), 
		definitionVars = list(), name = 'objective') {
		.Object@name <- name
		.Object@covariance <- covariance
		.Object@means <- means
		.Object@data <- data
		.Object@definitionVars <- definitionVars
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
		return(.Object)
})

setMethod("omxObjFunConvert", signature("MxMLObjective"), 
	function(.Object, flatModel, model) {
		name <- .Object@name
		covariance <- .Object@covariance
		means <- .Object@means
		data <- .Object@data
		covarianceIndex <- omxLocateIndex(flatModel, covariance, name)
		if(is.na(data)) {
			msg <- paste("The MxMLObjective", omxQuotes(name),
				"does not have a dataset associated with it in model",
				omxQuotes(flatModel@name))
			stop(msg, call.=FALSE)
		}
		if(!(length(means) == 1 && is.na(means))) {
			meansIndex <- omxLocateIndex(flatModel, means, name)
		} else {
			meansIndex <- as.integer(NA)
		}
		dIndex <- omxLocateIndex(flatModel, data, name)
		.Object@covariance <- covarianceIndex
		.Object@means <- meansIndex
		.Object@data <- dIndex
		.Object@definitionVars <- generateDefinitionList(flatModel)
		return(.Object)
})

mxMLObjective <- function(covariance, means = NA) {
	if (missing(covariance) || typeof(covariance) != "character") {
		stop("Covariance argument is not a string (the name of the expected covariance matrix)")
	}
	if (!(is.na(means) || typeof(means) == "character")) {
		stop("Means argument is not a string (the name of the expected means matrix)")
	}
	if (is.na(means)) means <- as.integer(NA)
	return(new("MxMLObjective", covariance, means))
}

displayMLObjective <- function(objective) {
	cat("MxMLObjective", omxQuotes(objective@name), '\n')
	cat("@covariance :", omxQuotes(objective@covariance), '\n')
	cat("@means :", omxQuotes(objective@means), '\n')
	if (length(objective@result) == 0) {
		cat("@result : empty\n")
	} else {
		cat("@result\n") 
		print(objective@result)
	}
	invisible(objective)
}


setMethod("print", "MxMLObjective", function(x,...) { 
	displayMLObjective(x) 
})

setMethod("show", "MxMLObjective", function(object) { 
	displayMLObjective(object) 
})
