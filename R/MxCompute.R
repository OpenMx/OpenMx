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
	   id = "integer",
	   "VIRTUAL"),
	 contains = "MxBaseNamed")

setClassUnion("MxCompute", c("NULL", "MxBaseCompute"))

setGeneric("convertForBackend",
	function(.Object, flatModel, model) {
		return(standardGeneric("convertForBackend"))
	})

setGeneric("assignId",
	function(.Object, id) {
		return(standardGeneric("assignId"))
	})

setMethod("assignId", signature("MxBaseCompute"),
	function(.Object, id) {
		.Object@id <- id
		.Object
	})

setGeneric("getFreeVarGroup",
	function(.Object) {
		return(standardGeneric("getFreeVarGroup"))
	})

setMethod("getFreeVarGroup", signature("MxBaseCompute"),
	function(.Object) {
		list()
	})

#----------------------------------------------------

setClass(Class = "MxComputeOperation",
	 contains = "MxBaseCompute",
	 representation = representation(
	   free.set = "MxOptionalChar"))

setMethod("qualifyNames", signature("MxComputeOperation"),
	function(.Object, modelname, namespace) {
		.Object@name <- imxIdentifier(modelname, .Object@name)
		.Object@free.set <- imxConvertIdentifier(.Object@free.set, modelname, namespace)
		.Object
	})

setMethod("getFreeVarGroup", signature("MxComputeOperation"),
	function(.Object) {
		if (length(.Object@free.set)) {
			list(.Object@id, .Object@free.set)
		} else {
			list()
		}
	})

setMethod("convertForBackend", signature("MxComputeOperation"),
	function(.Object, flatModel, model) {
		name <- .Object@name
		.Object
	})

#----------------------------------------------------

setClass(Class = "MxComputeOnce",
	 contains = "MxComputeOperation",
	 representation = representation(
	   what = "MxCharOrNumber",
	   context = "character",
	   gradient = "logical",
	   hessian = "logical",
	   start = "logical"))

setMethod("qualifyNames", signature("MxComputeOnce"),
	function(.Object, modelname, namespace) {
		.Object <- callNextMethod();
		for (sl in c('what')) {
			slot(.Object, sl) <- imxConvertIdentifier(slot(.Object, sl), modelname, namespace)
		}
		.Object
	})

setMethod("convertForBackend", signature("MxComputeOnce"),
	function(.Object, flatModel, model) {
		.Object <- callNextMethod();
		name <- .Object@name
		if (any(!is.integer(.Object@what))) {
			expNum <- match(.Object@what, names(flatModel@expectations))
			algNum <- match(.Object@what, append(names(flatModel@algebras),
							     names(flatModel@fitfunctions)))
			if (any(is.na(expNum)) && any(is.na(algNum))) {
				stop("Can only apply MxComputeOnce to MxAlgebra or MxExpectation")
			}
			if (!any(is.na(expNum))) {
					# Usually negative numbers indicate matrices; not here
				.Object@what <- - expNum
			} else {
				if (any(algNum > length(flatModel@algebras)) && length(algNum) > 1) {
					stop("MxComputeOnce cannot evaluate more than 1 fit function")
				}
				.Object@what <- algNum - 1L
			}
		}
		if (length(.Object@what) == 0) warning("MxComputeOnce with nothing will have no effect")
		if (.Object@start && .Object@what < 0) stop("start is only used with algebra or a fitfunction")
		.Object
	})

setMethod("initialize", "MxComputeOnce",
	  function(.Object, what, free.set, context, gradient, hessian, start) {
		  .Object@name <- 'compute'
		  .Object@what <- what
		  .Object@free.set <- free.set
		  .Object@context <- context
		  .Object@gradient <- gradient
		  .Object@hessian <- hessian
		  .Object@start <- start
		  .Object
	  })

mxComputeOnce <- function(what, free.set=NULL, context=character(0), gradient=FALSE,
			  hessian=FALSE, start=FALSE) {
	new("MxComputeOnce", what, free.set, context, gradient, hessian, start)
}

#----------------------------------------------------

setClass(Class = "MxComputeGradientDescent",
	 contains = "MxComputeOperation",
	 representation = representation(
	   start = "logical",
	   useGradient = "logical",
	   fitfunction = "MxCharOrNumber",
	   engine = "character",
	   verbose = "logical"))

setMethod("qualifyNames", signature("MxComputeGradientDescent"),
	function(.Object, modelname, namespace) {
		.Object <- callNextMethod();
		for (sl in c('fitfunction')) {
			slot(.Object, sl) <- imxConvertIdentifier(slot(.Object, sl), modelname, namespace)
		}
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
	  function(.Object, free.set, engine, fit, useGradient, start, verbose) {
		  .Object@name <- 'compute'
		  .Object@free.set <- free.set
		  .Object@fitfunction <- fit
		  .Object@start <- start
		  .Object@engine <- engine
		  .Object@useGradient <- useGradient
		  .Object@verbose <- verbose
		  .Object
	  })

mxComputeGradientDescent <- function(type=NULL, free.set=NULL, useGradient=as.logical(NA),
				     engine=NULL, fitfunction='fitfunction', start=FALSE, verbose=FALSE) {
# What to do with 'type'?
#	if (length(type) != 1) stop("Specific 1 compute type")

	if (is.null(engine)) engine <- as.character(NA)

	new("MxComputeGradientDescent", free.set, engine, fitfunction, useGradient, start, verbose)
}

#----------------------------------------------------

setClass(Class = "MxComputeNewtonRaphson",
	 contains = "MxComputeOperation",
	 representation = representation(
	   start = "logical",
	   fitfunction = "MxCharOrNumber",
	   maxIter = "integer",
	   tolerance = "numeric",
	   verbose = "logical"))

setMethod("qualifyNames", signature("MxComputeNewtonRaphson"),
	function(.Object, modelname, namespace) {
		.Object <- callNextMethod();
		for (sl in c('fitfunction')) {
			slot(.Object, sl) <- imxConvertIdentifier(slot(.Object, sl), modelname, namespace)
		}
		.Object
	})

setMethod("convertForBackend", signature("MxComputeNewtonRaphson"),
	function(.Object, flatModel, model) {
		.Object <- callNextMethod();
		name <- .Object@name
		if (is.character(.Object@fitfunction)) {
			.Object@fitfunction <- imxLocateIndex(flatModel, .Object@fitfunction, name)
		}
		.Object
	})

setMethod("initialize", "MxComputeNewtonRaphson",
	  function(.Object, free.set, fit, maxIter, tolerance, start, verbose) {
		  .Object@name <- 'compute'
		  .Object@free.set <- free.set
		  .Object@fitfunction <- fit
		  .Object@start <- start
		  .Object@maxIter <- maxIter
		  .Object@tolerance <- tolerance
		  .Object@verbose <- verbose
		  .Object
	  })

mxComputeNewtonRaphson <- function(type, free.set=NULL, start=FALSE,
				   fitfunction='fitfunction', maxIter = 100L, tolerance=1e-7,
				   verbose=FALSE) {

	new("MxComputeNewtonRaphson", free.set, fitfunction, maxIter, tolerance, start, verbose)
}

#----------------------------------------------------

setClass(Class = "MxComputeSteps",
	 contains = "MxBaseCompute",
	 representation = representation(
	   steps = "list"))

setMethod("getFreeVarGroup", signature("MxComputeSteps"),
	function(.Object) {
		result <- list()
		for (step in .Object@steps) {
			got <- getFreeVarGroup(step)
			if (length(got)) result <- append(result, got)
		}
		result
	})

setMethod("assignId", signature("MxComputeSteps"),
	function(.Object, id) {
		steps <- .Object@steps
		for (sx in 1:length(steps)) {
			steps[[sx]] <- assignId(steps[[sx]], id)
			id <- steps[[sx]]@id + 1L
		}
		.Object@steps <- steps
		.Object@id <- id
		.Object
	})

setMethod("qualifyNames", signature("MxComputeSteps"),
	function(.Object, modelname, namespace) {
		.Object@name <- imxIdentifier(modelname, .Object@name)
		.Object@steps <- lapply(.Object@steps, function (c) qualifyNames(c, modelname, namespace))
		.Object
	})

setMethod("convertForBackend", signature("MxComputeSteps"),
	function(.Object, flatModel, model) {
		.Object@steps <- lapply(.Object@steps, function (c) convertForBackend(c, flatModel, model))
		.Object
	})

#----------------------------------------------------

setClass(Class = "MxComputeIterate",
	 contains = "MxComputeSteps",
	 representation = representation(
	   maxIter = "integer",
	   tolerance = "numeric",
	   verbose = "logical"))

setMethod("initialize", "MxComputeIterate",
	  function(.Object, steps, maxIter, tolerance, verbose) {
		  .Object@name <- 'compute'
		  .Object@steps <- steps
		  .Object@maxIter <- maxIter
		  .Object@tolerance <- tolerance
		  .Object@verbose <- verbose
		  .Object
	  })

mxComputeIterate <- function(steps, maxIter=500L, tolerance=1e-4, verbose=FALSE) {
	new("MxComputeIterate", steps=steps, maxIter=maxIter, tolerance=tolerance, verbose)
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
		.Object <- callNextMethod();
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
	  function(.Object, free.set, fit, want.se) {
		  .Object@name <- 'compute'
		  .Object@free.set <- free.set
		  .Object@fitfunction <- fit
		  .Object@se <- want.se
		  .Object
	  })

mxComputeEstimatedHessian <- function(free.set=NULL, fitfunction='fitfunction', want.se=TRUE) {
	new("MxComputeEstimatedHessian", free.set, fitfunction, want.se)
}

#----------------------------------------------------

setClass(Class = "MxComputeSequence",
	 contains = "MxComputeSteps")

setMethod("initialize", "MxComputeSequence",
	  function(.Object, steps) {
		  .Object@name <- 'compute'
		  .Object@steps <- steps
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
	cat("@id :", opt@id, '\n')
	cat("@free.set :", omxQuotes(opt@free.set), '\n')
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
