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


setClass(Class = "MxFIMLObjective",
	representation = representation(
		covariance = "MxCharOrNumber",
		means = "MxCharOrNumber",
		definitionVars = "list"),
	contains = "MxBaseObjective")

setMethod("initialize", "MxFIMLObjective",
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

setMethod("omxObjFunNamespace", signature("MxFIMLObjective"), 
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


setMethod("omxObjFunConvert", signature("MxFIMLObjective"), 
	function(.Object, flatModel, model) {
		name <- .Object@name
		if(is.na(.Object@data)) {
			msg <- paste("The MxFIMLObjective", omxQuotes(name),
				"does not have a dataset associated with it in model",
				omxQuotes(flatModel@name))
			stop(msg, call.=FALSE)
		}
		if(flatModel@datasets[[.Object@data]]@type != 'raw') {
			msg <- paste("The dataset associated with MxFIMLObjective", 
				omxQuotes(name), "in model", 
				omxQuotes(flatModel@name), "is not raw data.")
			stop(msg, call.=FALSE)
		}
		if(!is.na(.Object@means)) {
			.Object@means <- omxLocateIndex(flatModel, .Object@means, name)
		}
		.Object@covariance <- omxLocateIndex(flatModel, .Object@covariance, name)
		.Object@data <- omxLocateIndex(flatModel, .Object@data, name)
		.Object@definitionVars <- generateDefinitionList(flatModel)
		return(.Object)
})

mxFIMLObjective <- function(covariance, means) {
	if (missing(covariance) || typeof(covariance) != "character") {
		stop("Covariance argument is not a string (the name of the expected covariance matrix)")
	}
	if (missing(means) || typeof(means) != "character") {
		stop("Means argument is not a string (the name of the expected means vector)")
	}
	return(new("MxFIMLObjective", covariance, means))
}

displayFIMLObjective <- function(objective) {
	cat("MxFIMLObjective", omxQuotes(objective@name), '\n')
	cat("covariance :", omxQuotes(objective@covariance), '\n')
	cat("means :", omxQuotes(objective@means), '\n')
	if (length(objective@result) == 0) {
		cat("Result : empty\n")
	} else {
		cat("Result : \n") 
		print(objective@result)
	}
	invisible(objective)
}


setMethod("print", "MxFIMLObjective", function(x,...) { 
	displayFIMLObjective(x) 
})

setMethod("show", "MxFIMLObjective", function(object) { 
	displayFIMLObjective(object) 
})