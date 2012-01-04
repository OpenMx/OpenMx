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


setClass(Class = "MxRObjective",
	representation = representation(
		objfun = "function",
		model = "MxModel",
		flatModel = "MxFlatModel",
		parameters = "list",
		state = "list"),
	contains = "MxBaseObjective")

setMethod("initialize", "MxRObjective",
	function(.Object, objfun, state, data = as.integer(NA), 
		name = 'objective') {
		.Object@objfun <- objfun
		.Object@name <- name
		.Object@data <- data
		.Object@state <- state
		return(.Object)
	}
)

setMethod("genericObjFunConvert", signature("MxRObjective"), 
	function(.Object, flatModel, model, defVars) {
		.Object@model <- model
		.Object@flatModel <- flatModel
		.Object@parameters <- generateParameterList(flatModel)
		return(.Object)
})

setMethod("genericObjFunNamespace", signature("MxRObjective"), 
	function(.Object, modelname, namespace) {
		.Object@name <- imxIdentifier(modelname, .Object@name)
		return(.Object)
})

setMethod("genericObjModelConvert", signature("MxRObjective"), 
	function(.Object, job, model, namespace, labelsData, flatJob) {
		job@.forcesequential <- TRUE
		job@.newobjects <- FALSE
		job@.newobjective <- FALSE
		job@.newtree <- FALSE
		return(list(job, flatJob))
})

mxRObjective <- function(objfun, ...) {
	if (!is.function(objfun)) {
		stop("First argument 'objfun' must be of type function")
	}
	if (length(formals(objfun)) != 2) {
		stop("The objective function must take exactly two arguments: a model and a persistant state")
	}
	state <- list(...)
	return(new("MxRObjective", objfun, state))
}

displayRObjective <- function(objective) {
	cat("MxRObjective", omxQuotes(objective@name), '\n')
	cat("@objfun (objective function) \n")
	print(objective@objfun)
	if (length(objective@result) == 0) {
		cat("@result: (not yet computed) ")
	} else {
		cat("@result:\n")
	}
	print(objective@result)
	invisible(objective)
}


setMethod("print", "MxRObjective", function(x,...) { 
	displayRObjective(x) 
})

setMethod("show", "MxRObjective", function(object) { 
	displayRObjective(object) 
})
