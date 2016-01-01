#
#   Copyright 2007-2016 The OpenMx Project
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


setClass(Class = "MxFitFunctionR",
	representation = representation(
		fitfun = "function",
	    units = "character",
		model = "MxModel",
		flatModel = "MxFlatModel",
		state = "list"),
	contains = "MxBaseFitFunction")

setMethod("initialize", "MxFitFunctionR",
	function(.Object, fitfun, state,  units, name = 'fitfunction') {
		.Object@fitfun <- fitfun
		.Object@name <- name
		.Object@state <- state
		.Object@expectation <- as.integer(NA)
		.Object@units <- units
		return(.Object)
	}
)

setMethod("genericFitFunConvert", signature("MxFitFunctionR"), 
	function(.Object, flatModel, model, labelsData, dependencies) {
		.Object@model <- model
		.Object@flatModel <- flatModel
		return(.Object)
})

setMethod("qualifyNames", signature("MxFitFunctionR"), 
	function(.Object, modelname, namespace) {
		.Object@name <- imxIdentifier(modelname, .Object@name)
		return(.Object)
})

mxFitFunctionR <- function(fitfun, ..., units="-2lnL") {
	if (!is.function(fitfun)) {
		stop("First argument 'fitfun' must be of type function")
	}
	if (length(formals(fitfun)) != 2) {
		stop("The fit function must take exactly two arguments: a model and a persistant state")
	}
	state <- list(...)
	return(new("MxFitFunctionR", fitfun, state, units))
}

displayRFitFun <- function(fitfunction) {
	cat("MxFitFunctionR", omxQuotes(fitfunction@name), '\n')
	cat("$fitfun (fitness function) \n")
	cat("$units: ", omxQuotes(fitfunction@units), '\n')
	print(fitfunction@fitfun)
	if (length(fitfunction@result) == 0) {
		cat("$result: (not yet computed) ")
	} else {
		cat("$result:\n")
	}
	print(fitfunction@result)
	invisible(fitfunction)
}


setMethod("print", "MxFitFunctionR", function(x,...) { 
	displayRFitFun(x) 
})

setMethod("show", "MxFitFunctionR", function(object) { 
	displayRFitFun(object) 
})
