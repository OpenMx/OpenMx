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
	   verbose = "integer",
	   context = "character",
	   maxAbsChange = "logical",
	   fit = "logical",
	   gradient = "logical",
	   hessian = "logical",
	     information = "logical",
	   ihessian = "logical",
	   hgprod = "logical"))

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
				stop(paste("Can only apply MxComputeOnce to MxAlgebra or MxExpectation not",
					   deparse(.Object@what)))
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
		if (all(.Object@what >= 0) && !.Object@maxAbsChange && !.Object@fit && !.Object@gradient &&
			    !.Object@hessian && !.Object@information && !.Object@ihessian && !.Object@hgprod) {
			warning("MxComputeOnce with no action")
		}
		.Object
	})

setMethod("initialize", "MxComputeOnce",
	  function(.Object, what, free.set, context, maxAbsChange, fit, gradient,
		   hessian, information, ihessian, hgprod, verbose) {
		  .Object@name <- 'compute'
		  .Object@what <- what
		  .Object@verbose = verbose
		  .Object@free.set <- free.set
		  .Object@context <- context
		  .Object@maxAbsChange <- maxAbsChange
		  .Object@fit <- fit
		  .Object@gradient <- gradient
		  .Object@hessian <- hessian
		  .Object@information <- information
		  .Object@ihessian <- ihessian
		  .Object@hgprod <- hgprod
		  .Object
	  })

mxComputeOnce <- function(what, free.set=NULL, context=character(0),
			  maxAbsChange=FALSE, fit=FALSE, gradient=FALSE,
			  hessian=FALSE, information=FALSE, ihessian=FALSE, hgprod=FALSE, verbose=0L) {
	new("MxComputeOnce", what, free.set, context, maxAbsChange, fit, gradient,
	    hessian, information, ihessian, hgprod, verbose)
}

#----------------------------------------------------

setClass(Class = "MxComputeGradientDescent",
	 contains = "MxComputeOperation",
	 representation = representation(
	   useGradient = "MxOptionalLogical",
	   fitfunction = "MxCharOrNumber",
	   engine = "character",
	   verbose = "integer"))

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
	  function(.Object, free.set, engine, fit, useGradient, verbose) {
		  .Object@name <- 'compute'
		  .Object@free.set <- free.set
		  .Object@fitfunction <- fit
		  .Object@engine <- engine
		  .Object@useGradient <- useGradient
		  .Object@verbose <- verbose
		  .Object
	  })

##' imxHasNPSOL
##'
##' @return
##' Returns TRUE if the NPSOL proprietary optimizer is compiled and
##' linked with OpenMx. Otherwise FALSE.
imxHasNPSOL <- function() .Call(hasNPSOL_wrapper)

mxComputeGradientDescent <- function(type=NULL, free.set=NULL, useGradient=NULL,
				     engine=NULL, fitfunction='fitfunction', verbose=0L) {
# What to do with 'type'?

	if (missing(engine)) {
		engine <- options()$mxOptions[["Default optimizer"]]
	}

	new("MxComputeGradientDescent", free.set, engine, fitfunction, useGradient, verbose)
}

#----------------------------------------------------

setClass(Class = "MxComputeNewtonRaphson",
	 contains = "MxComputeOperation",
	 representation = representation(
	   fitfunction = "MxCharOrNumber",
	   maxIter = "integer",
	   tolerance = "numeric",
	   verbose = "integer",
	   carefully = "logical"))

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
	  function(.Object, free.set, fit, maxIter, tolerance, verbose, carefully) {
		  .Object@name <- 'compute'
		  .Object@free.set <- free.set
		  .Object@fitfunction <- fit
		  .Object@maxIter <- maxIter
		  .Object@tolerance <- tolerance
		  .Object@verbose <- verbose
		  .Object@carefully <- carefully
		  .Object
	  })

mxComputeNewtonRaphson <- function(type, free.set=NULL,   # remove type arg TODO
				   fitfunction='fitfunction', maxIter = 100L, tolerance=1e-7,
				   verbose=0L, carefully=FALSE) {

	new("MxComputeNewtonRaphson", free.set, fitfunction, maxIter, tolerance, verbose, carefully)
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
	   verbose = "integer"))

setMethod("initialize", "MxComputeIterate",
	  function(.Object, steps, maxIter, tolerance, verbose) {
		  .Object@name <- 'compute'
		  .Object@steps <- steps
		  .Object@maxIter <- maxIter
		  .Object@tolerance <- tolerance
		  .Object@verbose <- verbose
		  .Object
	  })

mxComputeIterate <- function(steps, maxIter=500L, tolerance=1e-4, verbose=0L) {
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

setClass(Class = "MxComputeEM",
	 contains = "MxComputeOperation",
	 representation = representation(
	     what = "MxCharOrNumber",
	     mstep.fit = "MxCompute",
	     fit = "MxCompute",
	     maxIter = "integer",
	     tolerance = "numeric",
	     verbose = "integer",
	     ramsay="logical",
	     information="logical"))

setMethod("assignId", signature("MxComputeEM"),
	function(.Object, id) {
		.Object@mstep.fit <- assignId(.Object@mstep.fit, id)
		.Object@fit <- assignId(.Object@fit, id + 1L)
		.Object@id <- id + 2L
		.Object
	})

setMethod("getFreeVarGroup", signature("MxComputeEM"),
	function(.Object) {
		result <- callNextMethod();
		for (step in c(.Object@mstep.fit, .Object@fit)) {
			got <- getFreeVarGroup(step)
			if (length(got)) result <- append(result, got)
		}
		result
	})

setMethod("qualifyNames", signature("MxComputeEM"),
	function(.Object, modelname, namespace) {
		.Object <- callNextMethod();
		for (sl in c('what')) {
			slot(.Object, sl) <- imxConvertIdentifier(slot(.Object, sl), modelname, namespace)
		}
		for (sl in c('mstep.fit', 'fit')) {
			slot(.Object, sl) <- qualifyNames(slot(.Object, sl), modelname, namespace)
		}
		.Object
	})

setMethod("convertForBackend", signature("MxComputeEM"),
	function(.Object, flatModel, model) {
		.Object <- callNextMethod();
		name <- .Object@name
		if (any(!is.integer(.Object@what))) {
			expNum <- match(.Object@what, names(flatModel@expectations))
			if (any(is.na(expNum))) {
				stop(paste("Can only apply MxComputeEM to MxExpectation",
					   deparse(.Object@what)))
			}
			.Object@what <- expNum - 1L
		}
		if (length(.Object@what) == 0) warning("MxComputeEM with nothing will have no effect")
		for (sl in c('mstep.fit', 'fit')) {
			slot(.Object, sl) <- convertForBackend(slot(.Object, sl), flatModel, model)
		}
		.Object
	})

setMethod("initialize", "MxComputeEM",
	  function(.Object, what, mstep.fit, fit, maxIter, tolerance, verbose, ramsay, information) {
		  .Object@name <- 'compute'
		  .Object@what <- what
		  .Object@mstep.fit <- mstep.fit
		  .Object@fit <- fit
		  .Object@maxIter <- maxIter
		  .Object@tolerance <- tolerance
		  .Object@verbose <- verbose
		  .Object@ramsay <- ramsay
		  .Object@information <- information
		  .Object
	  })

mxComputeEM <- function(expectation, mstep.fit, fit, maxIter=500L, tolerance=1e-4,
			verbose=0L, ramsay=TRUE, information=FALSE) {
	new("MxComputeEM", what=expectation, mstep.fit, fit, maxIter=maxIter,
		tolerance=tolerance, verbose, ramsay, information)
}

displayMxComputeEM <- function(opt) {
	cat(class(opt), omxQuotes(opt@name), '\n')
	cat("@tolerance :", omxQuotes(opt@tolerance), '\n')
	cat("@maxIter :", omxQuotes(opt@maxIter), '\n')
	invisible(opt)
}

setMethod("print", "MxComputeEM", function(x, ...) displayMxComputeEM(x))
setMethod("show",  "MxComputeEM", function(object) displayMxComputeEM(object))

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
