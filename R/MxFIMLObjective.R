#
#   Copyright 2007-2009 The OpenMx Project
#
#   Licensed under the Apache License, Version 2.0 (the "License");
#   you may not use this file except in compliance with the License.
#   You may obtain a copy of the License at
# 
#        http://www.apache.org/licenses/LICENSE-2.0
# 
#    Unless required by applicable law or agreed to in writing, software
#    distributed under the License is distributed on an "AS IS" BASIS,
#    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#   See the License for the specific language governing permissions and
#   limitations under the License.


setClass(Class = "MxFIMLObjective",
	representation = representation(
		covariance = "MxCharOrNumber",
		means = "MxCharOrNumber",
		data = "numeric",
		definitionVars = "list"),
	contains = "MxBaseObjective")

setMethod("initialize", "MxFIMLObjective",
	function(.Object, name, covariance, means, 
				data = NA_real_, definitionVars = list()) {
		.Object@name <- name
		.Object@covariance <- covariance
		.Object@means <- means
		.Object@data <- data
		.Object@definitionVars <- definitionVars
		return(.Object)
	}
)

setMethod("omxObjFunConvert", signature("MxFIMLObjective", "MxFlatModel"), 
	function(.Object, model, definitions) {
		name <- .Object@name
		covariance <- .Object@covariance
		means <- .Object@means
		covarianceIndex <- omxLocateIndex(model, covariance, name)
		meansIndex <- omxLocateIndex(model, means, name)
		dIndex <- omxDataIndex(.Object@name, model@datasets)
		if (is.na(covarianceIndex)) {
			stop(paste("Could not find a matrix/algebra with name", 
			covariance, "in the model."))
		}
		if (is.na(meansIndex)) {
			stop(paste("Could not find a matrix/algebra with name", 
			means, "in the model."))
		}
		if (is.na(dIndex)) {
			stop(paste("Could not find a data set for objective", 
			.Object@name, "in the model."))
		}
		return(new("MxFIMLObjective", name, covarianceIndex, meansIndex, dIndex, definitions))
})

mxFIMLObjective <- function(covariance, means, name = omxUntitledName()) {
	if (typeof(name) != "character") {
		stop("Name argument is not a string (the name of the objective function)")
	}
	if (missing(covariance) || typeof(covariance) != "character") {
		stop("Covariance argument is not a string (the name of the expected covariance matrix)")
	}
	if (missing(means) || typeof(means) != "character") {
		stop("Means argument is not a string (the name of the expected means matrix)")
	}
	return(new("MxFIMLObjective", name, covariance, means))
}
