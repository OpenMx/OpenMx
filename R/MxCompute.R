#
#   Copyright 2013 The OpenMx Project
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

setClass(Class = "MxBaseCompute", 
	 representation = representation(
	   "VIRTUAL"),
	 contains = "MxBaseNamed")

setClassUnion("MxCompute", c("NULL", "MxBaseCompute"))

setGeneric("convertForBackend",
	function(.Object, flatModel, model) {
		return(standardGeneric("convertForBackend"))
	})

setMethod("qualifyNames", signature("MxBaseCompute"),
	function(.Object, modelname, namespace) {
		.Object@name <- imxIdentifier(modelname, .Object@name)
		.Object@fitfunction <- imxConvertIdentifier(.Object@fitfunction, modelname, namespace)
		.Object
	})

setMethod("convertForBackend", signature("MxBaseCompute"),
	function(.Object, flatModel, model) {
		name <- .Object@name
		if (is.character(.Object@fitfunction)) {
			.Object@fitfunction <- imxLocateIndex(flatModel, .Object@fitfunction, name)
		}
		fg <- match(.Object@free.group, flatModel@freeGroupNames)
		if (is.na(fg)) {
			stop(paste("Cannot find free group", .Object@free.group,
				   "in list of free groups:",
				   omxQuotes(flatModel@freeGroupNames)))
		} else {
			.Object@free.group <- fg - 1L
		}
		.Object
	})

setClass(Class = "MxComputeOperation",
	 contains = "MxBaseCompute",
	 representation = representation(
	   free.group = "MxCharOrNumber",
	   fitfunction = "MxCharOrNumber"))

setClass(Class = "MxComputeOnce",
	 contains = "MxComputeOperation")

setClass(Class = "MxComputeGradientDescent",
	 contains = "MxComputeOperation",
	 representation = representation(
	   type = "character",
	   engine = "character"))

setClass(Class = "MxComputeSequence",
	 contains = "MxBaseCompute",
	 representation = representation(
	   steps = "list"))

setClass(Class = "MxComputeEstimatedHessian",
	 contains = "MxComputeOperation",
	 representation = representation(
	   se = "logical"))

setMethod("initialize", "MxComputeOnce",
	  function(.Object, free.group, fit) {
		  .Object@name <- 'compute'
		  .Object@free.group <- free.group
		  .Object@fitfunction <- fit
		  .Object
	  })

mxComputeOnce <- function(free.group='default', fitfunction='fitfunction') {
	new("MxComputeOnce", free.group, fitfunction)
}

setMethod("initialize", "MxComputeSequence",
	  function(.Object, steps) {
		  .Object@name <- 'compute'
		  .Object@steps <- steps
		  .Object
	  })

setMethod("qualifyNames", signature("MxComputeSequence"),
	function(.Object, modelname, namespace) {
		.Object@name <- imxIdentifier(modelname, .Object@name)
		.Object@steps <- lapply(.Object@steps, function (c) qualifyNames(c, modelname, namespace))
		.Object
	})

setMethod("convertForBackend", signature("MxComputeSequence"),
	function(.Object, flatModel, model) {
		.Object@steps <- lapply(.Object@steps, function (c) convertForBackend(c, flatModel, model))
		.Object
	})

mxComputeSequence <- function(steps) {
	new("MxComputeSequence", steps=steps)
}

setMethod("initialize", "MxComputeEstimatedHessian",
	  function(.Object, free.group, fit, want.se) {
		  .Object@name <- 'compute'
		  .Object@free.group <- free.group
		  .Object@fitfunction <- fit
		  .Object@se <- want.se
		  .Object
	  })

mxComputeEstimatedHessian <- function(free.group='default', fitfunction='fitfunction', want.se=TRUE) {
	new("MxComputeEstimatedHessian", free.group, fitfunction, want.se)
}

setMethod("initialize", "MxComputeGradientDescent",
	  function(.Object, free.group, type, engine, fit) {
		  .Object@name <- 'compute'
		  .Object@free.group <- free.group
		  .Object@fitfunction <- fit
		  .Object@type <- type
		  .Object@engine <- engine
		  .Object
	  })

mxComputeGradientDescent <- function(type, free.group='default',
				     engine=NULL, fitfunction='fitfunction') {
#	if (length(type) != 1) stop("Specific 1 compute type")

	if (is.null(type)) type <- as.character(NA)
	if (is.null(engine)) engine <- as.character(NA)

	new("MxComputeGradientDescent", free.group, type, engine, fitfunction)
}

displayMxComputeSequence <- function(opt) {
	cat(class(opt), omxQuotes(opt@name), '\n')
	for (step in 1:length(opt@steps)) {
		cat("[[", step, "]] :", class(opt@steps[[step]]), '\n')
	}
	invisible(opt)
}

setMethod("print", "MxComputeSequence", function(x, ...) displayMxComputeSequence(x))
setMethod("show",  "MxComputeSequence", function(object) displayMxComputeSequence(object))

displayMxComputeOperation <- function(opt) {
	cat(class(opt), omxQuotes(opt@name), '\n')
	cat("@free.group :", omxQuotes(opt@free.group), '\n')
	cat("@fitfunction :", omxQuotes(opt@fitfunction), '\n')
	invisible(opt)
}

setMethod("print", "MxComputeOperation", function(x, ...) displayMxComputeOperation(x))
setMethod("show",  "MxComputeOperation", function(object) displayMxComputeOperation(object))

displayMxComputeGradientDescent <- function(opt) {
	cat("@type :", omxQuotes(opt@type), '\n')
	cat("@engine :", omxQuotes(opt@engine), '\n')
	invisible(opt)
}

setMethod("print", "MxComputeGradientDescent",
	  function(x, ...) { callNextMethod(); displayMxComputeGradientDescent(x) })
setMethod("show",  "MxComputeGradientDescent",
	  function(object) { callNextMethod(); displayMxComputeGradientDescent(object) })

convertComputes <- function(flatModel, model) {
	retval <- lapply(flatModel@computes, function(opt) {
		convertForBackend(opt, flatModel, model)
	})
	retval
}
