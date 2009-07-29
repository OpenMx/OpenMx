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


setClass(Class = "MxRObjective",
	representation = representation(
		objfun = "function",
		model = "MxModel",
		flatModel = "MxFlatModel",
		parameters = "list",
		env = "environment"),
	contains = "MxBaseObjective")

setMethod("initialize", "MxRObjective",
	function(.Object, objfun, data = as.numeric(NA), 
		name = 'objective') {
		.Object@objfun <- objfun
		.Object@name <- name
		.Object@data <- data
		return(.Object)
	}
)

setMethod("omxObjFunConvert", signature("MxRObjective"), 
	function(.Object, flatModel, model, env) {
		.Object@model <- model
		.Object@flatModel <- flatModel
		.Object@parameters <- generateParameterList(flatModel)
		.Object@env <- globalenv()
		return(.Object)
})

setMethod("omxObjFunNamespace", signature("MxRObjective"), 
	function(.Object, modelname, namespace) {
		.Object@name <- omxIdentifier(modelname, .Object@name)
		return(.Object)
})


mxRObjective <- function(objfun) {
	if (!is.function(objfun)) {
		stop("First argument 'objfun' must be of type function")
	}
	if (length(formals(objfun)) != 2) {
		stop("The objective function must take exactly two arguments: a model and a persistant state")
	}
	return(new("MxRObjective", objfun))
}

displayRObjective <- function(objective) {
	cat("MxRObjective", omxQuotes(objective@name), '\n')
	cat("Objective function : \n")
	print(objective@objfun)
	if (length(objective@result) == 0) {
		cat("Result : empty\n")
	} else {
		cat("Result : \n") 
		print(objective@result)
	}
	invisible(objective)
}


setMethod("print", "MxRObjective", function(x,...) { 
	displayRObjective(x) 
})

setMethod("show", "MxRObjective", function(object) { 
	displayRObjective(object) 
})
