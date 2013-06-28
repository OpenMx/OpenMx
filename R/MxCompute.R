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

#----------------------------------------------------

setClass(Class = "MxComputeOperation",
	 contains = "MxBaseCompute",
	 representation = representation(
	   free.group = "MxCharOrNumber"))

setMethod("qualifyNames", signature("MxComputeOperation"),
	function(.Object, modelname, namespace) {
		.Object@name <- imxIdentifier(modelname, .Object@name)
		.Object
	})

setMethod("convertForBackend", signature("MxComputeOperation"),
	function(.Object, flatModel, model) {
		name <- .Object@name
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

#----------------------------------------------------

setClass(Class = "MxComputeAssign",  # good name? or ComputeCopy?
	 contains = "MxComputeOperation",
	 representation = representation(
	   from = "MxCharOrNumber",
	   to = "MxCharOrNumber"))

setMethod("initialize", "MxComputeAssign",
	  function(.Object, from, to, free.group) {
		  .Object@name <- 'compute'
		  .Object@from <- from
		  .Object@to <- to
		  .Object@free.group <- free.group
		  .Object
	  })

setMethod("qualifyNames", signature("MxComputeAssign"),
	function(.Object, modelname, namespace) {
		.Object <- callNextMethod();
		.Object@from <- imxIdentifier(modelname, .Object@from)
		.Object@to <- imxIdentifier(modelname, .Object@to)
		.Object
	})

setMethod("convertForBackend", signature("MxComputeAssign"),
	function(.Object, flatModel, model) {
		.Object <- callNextMethod();
		name <- .Object@name
		for (sl in c('from', 'to')) {
			mat <- match(slot(.Object, sl), names(flatModel@matrices))
			if (any(is.na(mat))) {
				stop(paste("MxComputeAssign: cannot find",
					   omxQuotes(slot(.Object, sl)[is.na(mat)]),
					   "mentioned in slot '", sl, "'"))
			}
			slot(.Object, sl) <- -mat
		}
		.Object
	})

mxComputeAssign <- function(from, to, free.group="default") {
	if (length(from) != length(to)) {
		stop("Arguments 'from' and 'to' must be the same length")
	}
	new("MxComputeAssign", from=from, to=to, free.group=free.group)
}

#----------------------------------------------------

setClass(Class = "MxComputeOnce",
	 contains = "MxComputeOperation",
	 representation = representation(
	   fitfunction = "MxOptionalCharOrNumber",
	   expectation = "MxOptionalCharOrNumber",
	   context = "character"))

setMethod("qualifyNames", signature("MxComputeOnce"),
	function(.Object, modelname, namespace) {
		.Object@name <- imxIdentifier(modelname, .Object@name)
		.Object@fitfunction <- imxConvertIdentifier(.Object@fitfunction, modelname, namespace)
		.Object@expectation <- imxConvertIdentifier(.Object@expectation, modelname, namespace)
		.Object
	})

setMethod("convertForBackend", signature("MxComputeOnce"),
	function(.Object, flatModel, model) {
		.Object <- callNextMethod();
		name <- .Object@name
		if (is.character(.Object@fitfunction)) {
			.Object@fitfunction <- imxLocateIndex(flatModel, .Object@fitfunction, name)
		}
		if (is.character(.Object@expectation)) {
			.Object@expectation <- imxLocateIndex(flatModel, .Object@expectation, name)
		}
		.Object
	})

setMethod("initialize", "MxComputeOnce",
	  function(.Object, free.group, fit, expectation, context) {
		  .Object@name <- 'compute'
		  .Object@free.group <- free.group
		  if (!is.null(fit) && !is.null(expectation)) {
			  stop("Cannot evaluate a fitfunction and expectation simultaneously")
		  }
		  .Object@fitfunction <- fit
		  .Object@expectation <- expectation
		  .Object@context <- context
		  .Object
	  })

mxComputeOnce <- function(free.group='default', fitfunction=NULL, expectation=NULL, context=character(0)) {
	new("MxComputeOnce", free.group, fitfunction, expectation, context)
}

#----------------------------------------------------

setClass(Class = "MxComputeGradientDescent",
	 contains = "MxComputeOperation",
	 representation = representation(
	   fitfunction = "MxCharOrNumber",
	   engine = "character"))

setMethod("qualifyNames", signature("MxComputeGradientDescent"),
	function(.Object, modelname, namespace) {
		.Object@name <- imxIdentifier(modelname, .Object@name)
		.Object@fitfunction <- imxConvertIdentifier(.Object@fitfunction, modelname, namespace)
		.Object
	})

setMethod("convertForBackend", signature("MxComputeGradientDescent"),
	function(.Object, flatModel, model) {
		.Object <- callNextMethod();
		name <- .Object@name
		if (is.character(.Object@fitfunction)) {
			.Object@fitfunction <- imxLocateIndex(flatModel, .Object@fitfunction, name)
		}
		.Object
	})

setMethod("initialize", "MxComputeGradientDescent",
	  function(.Object, free.group, engine, fit) {
		  .Object@name <- 'compute'
		  .Object@free.group <- free.group
		  .Object@fitfunction <- fit
		  .Object@engine <- engine
		  .Object
	  })

mxComputeGradientDescent <- function(type, free.group='default',
				     engine=NULL, fitfunction='fitfunction') {
# What to do with 'type'?
#	if (length(type) != 1) stop("Specific 1 compute type")

	if (is.null(engine)) engine <- as.character(NA)

	new("MxComputeGradientDescent", free.group, engine, fitfunction)
}

#----------------------------------------------------

setClass(Class = "MxComputeIterate",
	 contains = "MxBaseCompute",
	 representation = representation(
	   steps = "list",
	   maxIter = "integer",
	   tolerance = "numeric"))

setMethod("initialize", "MxComputeIterate",
	  function(.Object, steps, maxIter, tolerance) {
		  .Object@name <- 'compute'
		  .Object@steps <- steps
		  .Object@maxIter <- maxIter
		  .Object@tolerance <- tolerance
		  .Object
	  })

setMethod("qualifyNames", signature("MxComputeIterate"),
	function(.Object, modelname, namespace) {
		.Object@name <- imxIdentifier(modelname, .Object@name)
		.Object@steps <- lapply(.Object@steps, function (c) qualifyNames(c, modelname, namespace))
		.Object
	})

setMethod("convertForBackend", signature("MxComputeIterate"),
	function(.Object, flatModel, model) {
		.Object@steps <- lapply(.Object@steps, function (c) convertForBackend(c, flatModel, model))
		.Object
	})

mxComputeIterate <- function(steps, maxIter=500L, tolerance=1e-4) {
	new("MxComputeIterate", steps=steps, maxIter=maxIter, tolerance=tolerance)
}

displayMxComputeIterate <- function(opt) {
	cat(class(opt), omxQuotes(opt@name), '\n')
	cat("@tolerance :", omxQuotes(opt@tolerance), '\n')
	cat("@maxIter :", omxQuotes(opt@maxIter), '\n')
	for (step in 1:length(opt@steps)) {
		cat("[[", step, "]] :", class(opt@steps[[step]]), '\n')
	}
	invisible(opt)
}

setMethod("print", "MxComputeIterate", function(x, ...) displayMxComputeIterate(x))
setMethod("show",  "MxComputeIterate", function(object) displayMxComputeIterate(object))

#----------------------------------------------------

setClass(Class = "MxComputeEstimatedHessian",
	 contains = "MxComputeOperation",
	 representation = representation(
	   fitfunction = "MxCharOrNumber",
	   se = "logical"))

setMethod("qualifyNames", signature("MxComputeEstimatedHessian"),
	function(.Object, modelname, namespace) {
		.Object@name <- imxIdentifier(modelname, .Object@name)
		.Object@fitfunction <- imxConvertIdentifier(.Object@fitfunction, modelname, namespace)
		.Object
	})

setMethod("convertForBackend", signature("MxComputeEstimatedHessian"),
	function(.Object, flatModel, model) {
		.Object <- callNextMethod();
		name <- .Object@name
		if (is.character(.Object@fitfunction)) {
			.Object@fitfunction <- imxLocateIndex(flatModel, .Object@fitfunction, name)
		}
		.Object
	})

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

#----------------------------------------------------

setClass(Class = "MxComputeSequence",
	 contains = "MxBaseCompute",
	 representation = representation(
	   steps = "list"))

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

displayMxComputeSequence <- function(opt) {
	cat(class(opt), omxQuotes(opt@name), '\n')
	for (step in 1:length(opt@steps)) {
		cat("[[", step, "]] :", class(opt@steps[[step]]), '\n')
	}
	invisible(opt)
}

setMethod("print", "MxComputeSequence", function(x, ...) displayMxComputeSequence(x))
setMethod("show",  "MxComputeSequence", function(object) displayMxComputeSequence(object))

#----------------------------------------------------

displayMxComputeOperation <- function(opt) {
	cat(class(opt), omxQuotes(opt@name), '\n')
	cat("@free.group :", omxQuotes(opt@free.group), '\n')
	invisible(opt)
}

setMethod("print", "MxComputeOperation", function(x, ...) displayMxComputeOperation(x))
setMethod("show",  "MxComputeOperation", function(object) displayMxComputeOperation(object))

displayMxComputeGradientDescent <- function(opt) {
	cat("@type :", omxQuotes(opt@type), '\n')
	cat("@engine :", omxQuotes(opt@engine), '\n')
	cat("@fitfunction :", omxQuotes(opt@fitfunction), '\n')
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
