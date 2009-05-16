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
		covariance <- .Object@covariance
		means <- .Object@means
		data <- .Object@data
		covarianceIndex <- omxLocateIndex(flatModel, covariance, name)
		if(is.na(data)) {
			msg <- paste("The MxFIMLObjective", omxQuotes(name),
				"does not have a dataset associated with it in model",
				omxQuotes(flatModel@name))
			stop(msg, call.=FALSE)
		}
		if(flatModel@datasets[[data]]@type != 'raw') {
			msg <- paste("The dataset associated with MxFIMLObjective", 
				omxQuotes(name), "in model", 
				omxQuotes(flatModel@name), "is not raw data.")
			stop(msg, call.=FALSE)
		}
		if(!is.na(means)) {
			meansIndex <- omxLocateIndex(flatModel, means, name)
		} else {
			meansIndex <- means
		}
		dIndex <- omxLocateIndex(flatModel, data, name)
		.Object@covariance <- covarianceIndex
		.Object@means <- meansIndex
		.Object@data <- dIndex
		.Object@definitionVars <- generateDefinitionList(flatModel)
		return(.Object)
})

mxFIMLObjective <- function(covariance, means = NA) {
	if (missing(covariance) || typeof(covariance) != "character") {
		stop("Covariance argument is not a string (the name of the expected covariance matrix)")
	}
	if (!(is.na(means) || typeof(means) == "character")) {
		stop("Means argument is not a string (the name of the expected means matrix)")
	}
	if (is.na(means)) means <- as.numeric(NA)
	return(new("MxFIMLObjective", covariance, means))
}
