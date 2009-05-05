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
		data <- .Object@data
		covarianceIndex <- omxLocateIndex(model, covariance, name)
		if(is.na(data)) {
			msg <- paste("The MxFIMLObjective", omxQuotes(name),
				"does not have a dataset associated with it in model",
				omxQuotes(model@name))
			stop(msg, call.=FALSE)
		}
		if(model@datasets[[data]]@type != 'raw') {
			msg <- paste("The dataset associated with MxFIMLObjective", 
				omxQuotes(name), "in model", 
				omxQuotes(model@name), "is not raw data.")
			stop(msg, call.=FALSE)
		}
		if(!is.na(means)) {
			meansIndex <- omxLocateIndex(model, means, name)
		} else {
			meansIndex <- means
		}
		dIndex <- omxLocateIndex(model, data, name)
		return(new("MxFIMLObjective", name, covarianceIndex, meansIndex, dIndex, definitions))
})

mxFIMLObjective <- function(covariance, means = NA, name = NA) {
	if (is.na(name)) name <- omxUntitledName()
	omxVerifyName(name)
	if (typeof(name) != "character") {
		stop("Name argument is not a string (the name of the objective function)")
	}
	if (missing(covariance) || typeof(covariance) != "character") {
		stop("Covariance argument is not a string (the name of the expected covariance matrix)")
	}
	if (!(is.na(means) || typeof(means) == "character")) {
		stop("Means argument is not a string (the name of the expected means matrix)")
	}
	if(is.na(means)) means <- NA_real_
	return(new("MxFIMLObjective", name, covariance, means))
}
