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

setClass(Class = "BaseCompute",
	 representation = representation(
	   id = "integer",
	     freeSet = "MxOptionalChar",
	     output = "list",
	     debug = "list",
	   "VIRTUAL"),
	 contains = "MxBaseNamed")

setClassUnion("MxCompute", c("NULL", "BaseCompute"))

setGeneric("displayCompute",
	   function(Ob, indent) {
		   return(standardGeneric("displayCompute"))
	   })

setMethod("displayCompute", signature("BaseCompute"),
	  function(Ob, indent) {
		  sp <- paste(rep('  ', indent), collapse="")
		  cat(sp, class(Ob), omxQuotes(Ob@name), '\n')
#		  cat(sp, "$id :", Ob@id, '\n')   # only of interest to developers and introduces visual noise
		  cat(sp, "$freeSet :", omxQuotes(Ob@freeSet), '\n')
		  if (length(Ob$output)) {
			  for (elem in names(Ob$output)) {
				  cat(sp, "$output[[", omxQuotes(elem), "]] :", Ob@output[[elem]], '\n')
			  }
		  }
		  if (length(Ob$debug)) {
			  cat(sp, "$debug :", Ob@debug, '\n')
			  for (elem in names(Ob$debug)) {
				  cat(sp, "$debug[[", omxQuotes(elem), "]] :", Ob@debug[[elem]], '\n')
			  }
		  }
		  invisible(Ob)
	  })

setMethod("print", "BaseCompute", function(x, ...) displayCompute(x, 1L))
setMethod("show",  "BaseCompute", function(object) displayCompute(object, 1L))

setGeneric("convertForBackend",
	function(.Object, flatModel, model) {
		return(standardGeneric("convertForBackend"))
	})

setMethod("convertForBackend", signature("BaseCompute"),
       function(.Object, flatModel, model) { .Object })

setGeneric("updateFromBackend",
	function(.Object, computes) {
		return(standardGeneric("updateFromBackend"))
	})

setGeneric("assignId",
	function(.Object, id, defaultFreeSet) {
		return(standardGeneric("assignId"))
	})

setMethod("assignId", signature("BaseCompute"),
	function(.Object, id, defaultFreeSet) {
		.Object@id <- id
		if (length(.Object@freeSet) == 1 && is.na(.Object@freeSet)) .Object@freeSet <- defaultFreeSet
		.Object
	})

setGeneric("getFreeVarGroup",
	function(.Object) {
		return(standardGeneric("getFreeVarGroup"))
	})

setMethod("getFreeVarGroup", signature("BaseCompute"),
	function(.Object) {
		if (length(.Object@freeSet) == 0 || (length(.Object@freeSet) == 1 && .Object@freeSet == '.')) {
			# none or all variables
		} else {
			list(.Object@id, .Object@freeSet)
		}
	})

setMethod("qualifyNames", signature("BaseCompute"),
	function(.Object, modelname, namespace) {
		.Object@name <- imxIdentifier(modelname, .Object@name)
		.Object@freeSet <- imxConvertIdentifier(.Object@freeSet, modelname, namespace)
		.Object
	})

setMethod("updateFromBackend", signature("BaseCompute"),
	function(.Object, computes) {
		if (length(computes)) {
			mystuff <- which(.Object@id == computes[seq(1,length(computes),2)])
			if (length(mystuff)) {
				got <- computes[[2 * mystuff]]
				for (sl in names(got)) {
					slot(.Object, sl) <- got[[sl]]
				}
			}
		}
		.Object
	})
    
setMethod("$", "BaseCompute", imxExtractSlot)

setReplaceMethod("$", "BaseCompute",
	function(x, name, value) {
		return(imxReplaceSlot(x, name, value, check=TRUE))
	}
)

#----------------------------------------------------

setClass(Class = "MxComputeOnce",
	 contains = "BaseCompute",
	 representation = representation(
	     from = "MxCharOrNumber",
	     what = "character",
	     how = "MxOptionalChar",
	     verbose = "integer",
	     .is.bestfit="logical"))

setMethod("qualifyNames", signature("MxComputeOnce"),
	function(.Object, modelname, namespace) {
		.Object <- callNextMethod()
		for (sl in c('from')) {
			slot(.Object, sl) <- imxConvertIdentifier(slot(.Object, sl), modelname, namespace)
		}
		.Object
	})

setMethod("convertForBackend", signature("MxComputeOnce"),
	function(.Object, flatModel, model) {
		name <- .Object@name
		if (any(!is.integer(.Object@from))) {
			expNum <- match(.Object@from, names(flatModel@expectations))
			algNum <- match(.Object@from, append(names(flatModel@algebras),
							     names(flatModel@fitfunctions)))
			if (any(is.na(expNum)) && any(is.na(algNum))) {
				stop(paste("Can only apply MxComputeOnce to MxAlgebra or MxExpectation not",
					   deparse(.Object@from)))
			}
			if (!any(is.na(expNum))) {
					# Usually negative numbers indicate matrices; not here
				.Object@from <- - expNum
			} else {
				if (any(algNum > length(flatModel@algebras)) && length(algNum) > 1) {
					stop("MxComputeOnce cannot evaluate more than 1 fit function")
				}
				.Object@from <- algNum - 1L
			}
		}
		.Object
	})

setMethod("initialize", "MxComputeOnce",
	  function(.Object, from, what, how, freeSet, verbose, .is.bestfit) {
		  .Object@name <- 'compute'
		  .Object@from <- from
		  .Object@what <- what
		  .Object@how <- how
		  .Object@freeSet <- freeSet
		  .Object@verbose = verbose
		  .Object@.is.bestfit <- .is.bestfit
		  .Object
	  })

##' Compute something once
##'
##' Some models are optimized for a sparse Hessian. Therefore, it can
##' be much more efficient to compute the inverse Hessian in
##' comparison to computing the Hessian and then inverting it.
##'
##' The information matrix is only valid when parameters are at the
##' maximum likelihood estimate. The information matrix is returned in
##' model$output$hessian. You cannot request both the information
##' matrix and the Hessian. The information matrix is invarient to the
##' sign of the log likelihood scale whereas the Hessian is not.
##' Use the \code{how} parameter to specify which approximation to use
##' (one of "default", "hessian", "sandwich", "bread", and "meat").
##'
##' @param from the object to perform the computation (a vector of expectation or algebra names)
##' @param what what to compute (default is "nothing")
##' @param how to compute it (optional)
##' @param ...  Not used.  Forces remaining arguments to be specified by name.
##' @param freeSet names of matrices containing free variables
##' @param verbose the level of debugging output
##' @param .is.bestfit do not use; for backward compatibility
##' @aliases
##' MxComputeOnce-class
##' @examples
##' data(demoOneFactor)
##' factorModel <- mxModel(name ="One Factor",
##'   mxMatrix(type="Full", nrow=5, ncol=1, free=TRUE, values=0.2, name="A"),
##'     mxMatrix(type="Symm", nrow=1, ncol=1, free=FALSE, values=1, name="L"),
##'     mxMatrix(type="Diag", nrow=5, ncol=5, free=TRUE, values=1, name="U"),
##'     mxAlgebra(expression=A %*% L %*% t(A) + U, name="R"),
##'     mxFitFunctionML(),mxExpectationNormal(covariance="R", dimnames=names(demoOneFactor)),
##'     mxData(observed=cov(demoOneFactor), type="cov", numObs=500),
##'     mxComputeOnce('fitfunction', 'fit'))
##' factorModelFit <- mxRun(factorModel)
##' factorModelFit$output$fit  # 972.15

mxComputeOnce <- function(from, what="nothing", how=NULL, ...,
			  freeSet=NA_character_, verbose=0L, .is.bestfit=FALSE) {
	garbageArguments <- list(...)
	if (length(garbageArguments) > 0) {
		stop("mxComputeOnce does not accept values for the '...' argument")
	}
	if (length(from) == 0) warning("mxComputeOnce from nothing will have no effect")
	new("MxComputeOnce", from, what, how, freeSet, verbose, .is.bestfit)
}

setMethod("displayCompute", signature("MxComputeOnce"),
	  function(Ob, indent) {
		  callNextMethod();
		  sp <- paste(rep('  ', indent), collapse="")
		  for (sl in c("from", "what", "how", "verbose")) {
			  slname <- paste("$", sl, sep="")
			  if (is.null(slot(Ob, sl))) next
			  if (is.character(slot(Ob, sl))) {
				  cat(sp, slname, ":", omxQuotes(slot(Ob, sl)), '\n')
			  } else {
				  cat(sp, slname, ":", slot(Ob, sl), '\n')
			  }
		  }
		  invisible(Ob)
	  })

#----------------------------------------------------

setClass(Class = "MxComputeGradientDescent",
	 contains = "BaseCompute",
	 representation = representation(
	   useGradient = "MxOptionalLogical",
	   fitfunction = "MxCharOrNumber",
	   engine = "character",
	     tolerance = "numeric",
	   verbose = "integer"))

setMethod("qualifyNames", signature("MxComputeGradientDescent"),
	function(.Object, modelname, namespace) {
		.Object <- callNextMethod()
		for (sl in c('fitfunction')) {
			slot(.Object, sl) <- imxConvertIdentifier(slot(.Object, sl), modelname, namespace)
		}
		.Object
	})

setMethod("convertForBackend", signature("MxComputeGradientDescent"),
	function(.Object, flatModel, model) {
		name <- .Object@name
		if (is.character(.Object@fitfunction)) {
			.Object@fitfunction <- imxLocateIndex(flatModel, .Object@fitfunction, .Object)
		}
		.Object
	})

setMethod("initialize", "MxComputeGradientDescent",
	  function(.Object, freeSet, engine, fit, useGradient, verbose, tolerance) {
		  .Object@name <- 'compute'
		  .Object@freeSet <- freeSet
		  .Object@fitfunction <- fit
		  .Object@engine <- engine
		  .Object@useGradient <- useGradient
		  .Object@verbose <- verbose
		  .Object@tolerance <- tolerance
		  .Object
	  })

##' imxHasNPSOL
##'
##' @return
##' Returns TRUE if the NPSOL proprietary optimizer is compiled and
##' linked with OpenMx. Otherwise FALSE.
imxHasNPSOL <- function() .Call(hasNPSOL_wrapper)

##' Optimize parameters using a gradient descent optimizer
##'
##' This optimizer does not require analytic derivatives of the fit
##' function. The open-source version of OpenMx only offers 1 choice,
##' CSOLNP (based on Ye, 1988).  The proprietary version of OpenMx
##' offers the choice of two optimizers, CSOLNP and NPSOL.
##'
##' @param freeSet names of matrices containing free variables
##' @param ...  Not used.  Forces remaining arguments to be specified by name.
##' @param engine specific NPSOL or CSOLNP
##' @param fitfunction name of the fitfunction (defaults to 'fitfunction')
##' @param verbose level of debugging output
##' @param tolerance how close to the optimum is close enough (also known as the optimality tolerance)
##' @param useGradient whether to use the analytic gradient (if available)
##' @aliases
##' MxComputeGradientDescent-class
##' @references Ye, Y. (1988). \emph{Interior algorithms for linear,
##' quadratic, and linearly constrained convex programming.}
##' (Unpublished doctoral dissertation.) Stanford University, CA.
##' @examples
##' data(demoOneFactor)
##' factorModel <- mxModel(name ="One Factor",
##'   mxMatrix(type="Full", nrow=5, ncol=1, free=FALSE, values=0.2, name="A"),
##'     mxMatrix(type="Symm", nrow=1, ncol=1, free=FALSE, values=1, name="L"),
##'     mxMatrix(type="Diag", nrow=5, ncol=5, free=TRUE, values=1, name="U"),
##'     mxAlgebra(expression=A %*% L %*% t(A) + U, name="R"),
##'   mxExpectationNormal(covariance="R", dimnames=names(demoOneFactor)),
##'   mxFitFunctionML(),
##'     mxData(observed=cov(demoOneFactor), type="cov", numObs=500),
##'      mxComputeSequence(steps=list(
##'      mxComputeGradientDescent(),
##'      mxComputeNumericDeriv(),
##'      mxComputeStandardError(),
##'      mxComputeHessianQuality()
##'     )))
##' factorModelFit <- mxRun(factorModel)
##' factorModelFit$output$conditionNumber # 29.5

mxComputeGradientDescent <- function(freeSet=NA_character_, ...,
				     engine=NULL, fitfunction='fitfunction', verbose=0L,
				     tolerance=NA_real_, useGradient=NULL) {

	garbageArguments <- list(...)
	if (length(garbageArguments) > 0) {
		stop("mxComputeGradientDescent does not accept values for the '...' argument")
	}
	if (missing(engine)) {
		engine <- options()$mxOptions[["Default optimizer"]]
	}

	new("MxComputeGradientDescent", freeSet, engine, fitfunction, useGradient, verbose,
	    tolerance)
}

setMethod("displayCompute", signature("MxComputeGradientDescent"),
	  function(Ob, indent) {
		  callNextMethod();
		  sp <- paste(rep('  ', indent), collapse="")
		  cat(sp, "$engine :", omxQuotes(Ob@engine), '\n')
		  cat(sp, "$fitfunction :", omxQuotes(Ob@fitfunction), '\n')
		  cat(sp, "$verbose :", Ob@verbose, '\n')
		  if (!is.na(Ob@tolerance)) {
			  cat(sp, "$tolerance :", Ob@tolerance, '\n')
		  }
		  if (!is.null(Ob@useGradient)) {
			  cat(sp, "$useGradient :", Ob@useGradient, '\n')
		  }
		  invisible(Ob)
	  })

#----------------------------------------------------

setClass(Class = "MxComputeConfidenceInterval",
	 contains = "BaseCompute",
	 representation = representation(
	     fitfunction = "MxCharOrNumber",
	     engine = "character",
	     tolerance = "numeric",
	     verbose = "integer"))

setMethod("qualifyNames", signature("MxComputeConfidenceInterval"),
	function(.Object, modelname, namespace) {
		.Object <- callNextMethod()
		for (sl in c('fitfunction')) {
			slot(.Object, sl) <- imxConvertIdentifier(slot(.Object, sl), modelname, namespace)
		}
		.Object
	})

setMethod("convertForBackend", signature("MxComputeConfidenceInterval"),
	function(.Object, flatModel, model) {
		name <- .Object@name
		if (is.character(.Object@fitfunction)) {
			.Object@fitfunction <- imxLocateIndex(flatModel, .Object@fitfunction, .Object)
		}
		.Object
	})

setMethod("initialize", "MxComputeConfidenceInterval",
	  function(.Object, freeSet, engine, fit, verbose, tolerance) {
		  .Object@name <- 'compute'
		  .Object@freeSet <- freeSet
		  .Object@fitfunction <- fit
		  .Object@engine <- engine
		  .Object@verbose <- verbose
		  .Object@tolerance <- tolerance
		  .Object
	  })

##' Find likelihood-based confidence intervals
##'
##' Add some description TODO
##'
##' @param freeSet names of matrices containing free variables
##' @param ...  Not used.  Forces remaining arguments to be specified by name.
##' @param engine specific NPSOL or CSOLNP
##' @param fitfunction name of the fitfunction (defaults to 'fitfunction')
##' @param verbose level of debugging output
##' @param tolerance how close to the optimum is close enough (also known as the optimality tolerance)
##' @aliases
##' MxComputeConfidenceInterval-class

mxComputeConfidenceInterval <- function(freeSet=NA_character_, ...,
					engine=NULL, fitfunction='fitfunction', verbose=0L,
					tolerance=NA_real_) {

	garbageArguments <- list(...)
	if (length(garbageArguments) > 0) {
		stop("mxComputeConfidenceInterval does not accept values for the '...' argument")
	}
	if (missing(engine)) {
		engine <- options()$mxOptions[["Default optimizer"]]
	}

	new("MxComputeConfidenceInterval", freeSet, engine, fitfunction, verbose, tolerance)
}

setMethod("displayCompute", signature("MxComputeConfidenceInterval"),
	  function(Ob, indent) {
		  callNextMethod();
		  sp <- paste(rep('  ', indent), collapse="")
		  for (sl in c("fitfunction", "engine", "tolerance", "verbose")) {
			  if (is.na(slot(Ob, sl))) next
			  slname <- paste("$", sl, sep="")
			  if (is.character(slot(Ob, sl))) {
				  cat(sp, slname, ":", omxQuotes(slot(Ob, sl)), '\n')
			  } else {
				  cat(sp, slname, ":", slot(Ob, sl), '\n')
			  }
		  }
		  invisible(Ob)
	  })

#----------------------------------------------------

setClass(Class = "MxComputeNewtonRaphson",
	 contains = "BaseCompute",
	 representation = representation(
	   fitfunction = "MxCharOrNumber",
	   maxIter = "integer",
	   tolerance = "numeric",
	   verbose = "integer"))

setMethod("qualifyNames", signature("MxComputeNewtonRaphson"),
	function(.Object, modelname, namespace) {
		.Object <- callNextMethod()
		for (sl in c('fitfunction')) {
			slot(.Object, sl) <- imxConvertIdentifier(slot(.Object, sl), modelname, namespace)
		}
		.Object
	})

setMethod("convertForBackend", signature("MxComputeNewtonRaphson"),
	function(.Object, flatModel, model) {
		name <- .Object@name
		if (is.character(.Object@fitfunction)) {
			.Object@fitfunction <- imxLocateIndex(flatModel, .Object@fitfunction, .Object)
		}
		.Object
	})

setMethod("initialize", "MxComputeNewtonRaphson",
	  function(.Object, freeSet, fit, maxIter, tolerance, verbose) {
		  .Object@name <- 'compute'
		  .Object@freeSet <- freeSet
		  .Object@fitfunction <- fit
		  .Object@maxIter <- maxIter
		  .Object@tolerance <- tolerance
		  .Object@verbose <- verbose
		  .Object
	  })

##' Optimize parameters using the Newton-Raphson algorithm
##'
##' This optimizer requires analytic 1st and 2nd derivatives of the
##' fit function. Ramsay (1975) is used to speed convergence. Ramsay
##' can be differentially applied to different groups of parameters.
##' Comprehensive diagnostics are available by increasing the verbose
##' level.
##'
##' @param freeSet names of matrices containing free variables
##' @param ...  Not used.  Forces remaining arguments to be specified by name.
##' @param fitfunction name of the fitfunction (defaults to 'fitfunction')
##' @param maxIter maximum number of iterations
##' @param tolerance optimization is considered converged when maximum absolute change in parameters is less than tolerance
##' @param verbose level of debugging output
##' @aliases
##' MxComputeNewtonRaphson-class
##' @references
##' Luenberger, D. G. & Ye, Y. (2008). \emph{Linear and nonlinear programming.} Springer.
##' 
##' Ramsay, J. O. (1975). Solving implicit equations in psychometric data analysis.
##' \emph{Psychometrika, 40}(3), 337-360.

mxComputeNewtonRaphson <- function(freeSet=NA_character_, ..., fitfunction='fitfunction', maxIter = 100L,
				   tolerance=1e-7, verbose=0L)
{
	garbageArguments <- list(...)
	if (length(garbageArguments) > 0) {
		stop("mxComputeNewtonRaphson does not accept values for the '...' argument")
	}

	new("MxComputeNewtonRaphson", freeSet, fitfunction, maxIter, tolerance, verbose)
}

setMethod("displayCompute", signature("MxComputeEM"),
	  function(Ob, indent) {
		  callNextMethod();
		  sp <- paste(rep('  ', indent), collapse="")
		  for (sl in c("fitfunction", "maxIter", "tolerance", "verbose")) {
			  slname <- paste("$", sl, sep="")
			  if (is.character(slot(Ob, sl))) {
				  cat(sp, slname, ":", omxQuotes(slot(Ob, sl)), '\n')
			  } else {
				  cat(sp, slname, ":", slot(Ob, sl), '\n')
			  }
		  }
		  invisible(Ob)
	  })

#----------------------------------------------------

setClass(Class = "ComputeSteps",
	 contains = "BaseCompute",
	 representation = representation(
	   steps = "list"))

setMethod("getFreeVarGroup", signature("ComputeSteps"),
	function(.Object) {
		result <- callNextMethod()
		for (step in .Object@steps) {
			got <- getFreeVarGroup(step)
			if (length(got)) result <- append(result, got)
		}
		result
	})

setMethod("assignId", signature("ComputeSteps"),
	function(.Object, id, defaultFreeSet) {
		if (length(.Object@freeSet) == 1 && is.na(.Object@freeSet)) .Object@freeSet <- defaultFreeSet
		defaultFreeSet <- .Object@freeSet
		steps <- .Object@steps
		if (length(steps)) for (sx in 1:length(steps)) {
			steps[[sx]] <- assignId(steps[[sx]], id, defaultFreeSet)
			id <- steps[[sx]]@id + 1L
		}
		.Object@steps <- steps
		.Object@id <- id
		.Object
	})

setMethod("qualifyNames", signature("ComputeSteps"),
	function(.Object, modelname, namespace) {
		.Object <- callNextMethod()
		.Object@name <- imxIdentifier(modelname, .Object@name)
		.Object@steps <- lapply(.Object@steps, function (c) qualifyNames(c, modelname, namespace))
		.Object
	})

setMethod("convertForBackend", signature("ComputeSteps"),
	function(.Object, flatModel, model) {
		.Object@steps <- lapply(.Object@steps, function (c) convertForBackend(c, flatModel, model))
		.Object
	})

setMethod("updateFromBackend", signature("ComputeSteps"),
	function(.Object, computes) {
		.Object <- callNextMethod()
		.Object@steps <- lapply(.Object@steps, function (c) updateFromBackend(c, computes))
		.Object
	})

#----------------------------------------------------

setClass(Class = "MxComputeIterate",
	 contains = "ComputeSteps",
	 representation = representation(
	   maxIter = "integer",
	   tolerance = "numeric",
	   verbose = "integer"))

setMethod("initialize", "MxComputeIterate",
	  function(.Object, steps, maxIter, tolerance, verbose, freeSet) {
		  .Object@name <- 'compute'
		  .Object@steps <- steps
		  .Object@maxIter <- maxIter
		  .Object@tolerance <- tolerance
		  .Object@verbose <- verbose
		  .Object@freeSet <- freeSet
		  .Object
	  })

##' Repeatedly invoke a series of compute objects until change is less than tolerance
##'
##' One step (typically the last) must compute the fit or maxAbsChange.
##'
##' @param steps a list of compute objects
##' @param ...  Not used.  Forces remaining arguments to be specified by name.
##' @param maxIter the maximum number of iterations
##' @param tolerance iterates until change is less than tolerance
##' @param verbose level of debugging output
##' @param freeSet Names of matrices containing free variables.
##' @aliases
##' MxComputeIterate-class
mxComputeIterate <- function(steps, ..., maxIter=500L, tolerance=1e-4, verbose=0L, freeSet=NA_character_) {
	garbageArguments <- list(...)
	if (length(garbageArguments) > 0) {
		stop("mxComputeIterate does not accept values for the '...' argument")
	}

	new("MxComputeIterate", steps=steps, maxIter=maxIter, tolerance=tolerance, verbose, freeSet)
}

setMethod("displayCompute", signature("MxComputeIterate"),
	  function(Ob, indent) {
		  callNextMethod();
		  sp <- paste(rep('  ', indent), collapse="")
		  cat(sp, "maxIter :", Ob@maxIter, '\n')
		  cat(sp, "tolerance :", Ob@tolerance, '\n')
		  cat(sp, "verbose :", Ob@verbose, '\n')
		  for (step in 1:length(Ob@steps)) {
			  cat(sp, "steps[[", step, "]] :", '\n')
			  displayCompute(Ob@steps[[step]], indent+1L)
		  }
		  invisible(Ob)
	  })

#----------------------------------------------------

setClass(Class = "MxComputeEM",
	 contains = "BaseCompute",
	 representation = representation(
	     expectation = "MxCharOrNumber",
	     predict = "character",
	     mstep = "MxCompute",
	     observedFit = "MxCharOrNumber",
	     maxIter = "integer",
	     tolerance = "numeric",
	     verbose = "integer",
	     accel="character",
	     information="character",
	     infoArgs="list"))

setMethod("assignId", signature("MxComputeEM"),
	function(.Object, id, defaultFreeSet) {
		if (length(.Object@freeSet) == 1 && is.na(.Object@freeSet)) .Object@freeSet <- defaultFreeSet
		defaultFreeSet <- .Object@freeSet
		for (sl in c('mstep')) {
			slot(.Object, sl) <- assignId(slot(.Object, sl), id, defaultFreeSet)
			id <- slot(.Object, sl)@id + 1L
		}
		.Object@id <- id 
		.Object
	})

setMethod("getFreeVarGroup", signature("MxComputeEM"),
	function(.Object) {
		result <- callNextMethod()
		for (step in c(.Object@mstep)) {
			got <- getFreeVarGroup(step)
			if (length(got)) result <- append(result, got)
		}
		result
	})

setMethod("qualifyNames", signature("MxComputeEM"),
	function(.Object, modelname, namespace) {
		.Object <- callNextMethod()
		for (sl in c('expectation', 'observedFit')) {
			slot(.Object, sl) <- imxConvertIdentifier(slot(.Object, sl), modelname, namespace)
		}
		for (sl in c('mstep')) {
			slot(.Object, sl) <- qualifyNames(slot(.Object, sl), modelname, namespace)
		}
		.Object@infoArgs$fitfunction <-
		    imxConvertIdentifier(.Object@infoArgs$fitfunction, modelname, namespace)
		.Object
	})

setMethod("convertForBackend", signature("MxComputeEM"),
	function(.Object, flatModel, model) {
		name <- .Object@name
		if (any(!is.integer(.Object@expectation))) {
			expNum <- match(.Object@expectation, names(flatModel@expectations))
			if (any(is.na(expNum))) {
				stop(paste("MxComputeEM: MxExpectation",
					   deparse(.Object@expectation), "not found"))
			}
			.Object@expectation <- expNum - 1L
		}
		if (length(.Object@expectation) == 0) warning("MxComputeEM with nothing will have no effect")
		if (length(.Object@observedFit) != 1) stop("MxComputeEM requires a single observedFit function")
		if (any(!is.integer(.Object@observedFit))) {
			algNum <- match(.Object@observedFit, append(names(flatModel@algebras),
								    names(flatModel@fitfunctions)))
			if (any(is.na(algNum))) {
				stop(paste("MxComputeEM: observedFit fit function", .Object@observedFit, "not found"))
			}
			.Object@observedFit <- algNum - 1L
		}
		for (sl in c('mstep')) {
			slot(.Object, sl) <- convertForBackend(slot(.Object, sl), flatModel, model)
		}
		fit <- match(.Object@infoArgs$fitfunction,
			     append(names(flatModel@algebras), names(flatModel@fitfunctions)))
		if (any(is.na(fit))) {
			stop(paste("ComputeEM: cannot find fitfunction",
				   omxQuotes(.Object@infoArgs$fitfunction[is.na(fit)]), "in infoArgs"))
		}
		.Object@infoArgs$fitfunction <-fit - 1L

		.Object
	})

setMethod("updateFromBackend", signature("MxComputeEM"),
	function(.Object, computes) {
		.Object <- callNextMethod()
		for (sl in c('mstep')) {
			slot(.Object, sl) <- updateFromBackend(slot(.Object, sl), computes)
		}
		.Object
	})

setMethod("initialize", "MxComputeEM",
	  function(.Object, expectation, predict, mstep, observedFit, maxIter, tolerance,
		   verbose, accel, information, freeSet, infoArgs) {
		  .Object@name <- 'compute'
		  .Object@expectation <- expectation
		  .Object@predict <- predict
		  .Object@mstep <- mstep
		  .Object@observedFit <- observedFit
		  .Object@maxIter <- maxIter
		  .Object@tolerance <- tolerance
		  .Object@verbose <- verbose
		  .Object@accel <- accel
		  .Object@information <- information
		  .Object@freeSet <- freeSet
		  .Object@infoArgs <- infoArgs
		  .Object
	  })

##' Fit a model using DLR's (1977) Expectation-Maximization (EM) algorithm
##'
##' The EM algorithm constitutes the following steps: Start with an
##' initial parameter vector. Predict the missing data to form a
##' completed data model. Optimize the completed data model to obtain
##' a new parameter vector. Repeat these steps until convergence
##' criteria are met.
##'
##' This compute plan does not work with any and all expectations. It
##' requires a special kind of expectation that can predict its
##' missing data to create a completed data model.
##'
##' The EM algorithm does not produce a parameter covariance matrixn
##' for standard errors. S-EM, an implementation of Meng & Rubin
##' (1991), is included.
##'
##' @param expectation a vector of expectation names
##' @param predict what to predict from the observed data (available options depend on the expectation)
##' @param mstep a compute plan to optimize the completed data model
##' @param observedFit the name of the observed data fit function (defaults to "fitfunction")
##' @param ...  Not used.  Forces remaining arguments to be specified by name.
##' @param maxIter maximum number of iterations
##' @param tolerance maximum change in fit to judge convergence
##' @param verbose level of diagnostic output
##' @param freeSet names of matrices containing free variables
##' @param accel name of acceleration method (defaults to "ramsay1975")
##' @param information name of information matrix approximation method
##' @param infoArgs arguments to control the information matrix method
##' @aliases
##' MxComputeEM-class
##' @references
##' Dempster, A. P., Laird, N. M., & Rubin, D. B. (1977). Maximum likelihood from
##' incomplete data via the EM algorithm. \emph{Journal of the Royal Statistical Society.
##' Series B (Methodological)}, 1–38.
##'
##' Meng, X.-L. & Rubin, D. B. (1991). Using EM to obtain asymptotic variance-covariance
##' matrices: The SEM algorithm. \emph{Journal of the American Statistical Association,
##' 86} (416), 899–909.
##' 
##' Ramsay, J. O. (1975). Solving implicit equations in psychometric data analysis.
##' \emph{Psychometrika, 40} (3), 337–360.
mxComputeEM <- function(expectation, predict, mstep, observedFit="fitfunction", ...,
			maxIter=500L, tolerance=1e-4, verbose=0L, freeSet=NA_character_,
			accel="ramsay1975", information=NA_character_, infoArgs=list()) {
	garbageArguments <- list(...)
	if (length(garbageArguments) > 0) {
		stop("mxComputeEM does not accept values for the '...' argument")
	}
	new("MxComputeEM", expectation, predict, mstep, observedFit, maxIter=maxIter,
	    tolerance=tolerance, verbose, accel, information, freeSet, infoArgs)
}

setMethod("displayCompute", signature("MxComputeEM"),
	  function(Ob, indent) {
		  callNextMethod();
		  sp <- paste(rep('  ', indent), collapse="")
		  cat(sp, "$expectation :", omxQuotes(Ob@expectation), '\n')
		  cat(sp, "$predict :", omxQuotes(Ob@predict), '\n')
		  cat(sp, "$mstep :", '\n')
		  displayCompute(Ob@mstep, indent+1L)
		  for (sl in c("observedFit", "maxIter", "tolerance", "verbose", "accel")) {
			  slname <- paste("$", sl, sep="")
			  if (is.character(slot(Ob, sl))) {
				  cat(sp, slname, ":", omxQuotes(slot(Ob, sl)), '\n')
			  } else {
				  cat(sp, slname, ":", slot(Ob, sl), '\n')
			  }
		  }
		  if (!is.na(Ob@information)) {
			  cat(sp, "$information :", Ob@information, '\n')
		  }
		  invisible(Ob)
	  })

#----------------------------------------------------

setClass(Class = "MxComputeNumericDeriv",
	 contains = "BaseCompute",
	 representation = representation(
	   fitfunction = "MxCharOrNumber",
	     parallel = "logical",
	     stepSize = "numeric",
	     iterations = "integer",
	     verbose="integer"))

setMethod("qualifyNames", signature("MxComputeNumericDeriv"),
	function(.Object, modelname, namespace) {
		.Object <- callNextMethod();
		.Object@fitfunction <- imxConvertIdentifier(.Object@fitfunction, modelname, namespace)
		.Object
	})

setMethod("convertForBackend", signature("MxComputeNumericDeriv"),
	function(.Object, flatModel, model) {
		name <- .Object@name
		if (is.character(.Object@fitfunction)) {
			.Object@fitfunction <- imxLocateIndex(flatModel, .Object@fitfunction, .Object)
		}
		.Object
	})

setMethod("initialize", "MxComputeNumericDeriv",
	  function(.Object, freeSet, fit, parallel, stepSize, iterations, verbose) {
		  .Object@name <- 'compute'
		  .Object@freeSet <- freeSet
		  .Object@fitfunction <- fit
		  .Object@parallel <- parallel
		  .Object@stepSize <- stepSize
		  .Object@iterations <- iterations
		  .Object@verbose <- verbose
		  .Object
	  })

##' Numerically estimate Hessian using Richardson extrapolation
##'
##' For N free parameters, Richardson extrapolation requires
##' (iterations * (N^2 + N)) function evaluations.
##' 
##' The implementation is closely based on the numDeriv R package.
##' 
##' @param freeSet names of matrices containing free variables
##' @param ...  Not used.  Forces remaining arguments to be specified by name.
##' @param fitfunction name of the fitfunction (defaults to 'fitfunction')
##' @param parallel whether to evaluate the fitfunction in parallel (defaults to TRUE)
##' @param stepSize starting set size (defaults to 0.0001)
##' @param iterations number of Richardson extrapolation iterations (defaults to 4L)
##' @param verbose Level of debugging output.
##' @aliases
##' MxComputeNumericDeriv-class
##' @examples
##' library(OpenMx)
##' data(demoOneFactor)
##' factorModel <- mxModel(name ="One Factor",
##' 	mxMatrix(type = "Full", nrow = 5, ncol = 1, free = FALSE, values = .2, name = "A"), 
##' 	mxMatrix(type = "Symm", nrow = 1, ncol = 1, free = FALSE, values = 1 , name = "L"), 
##' 	mxMatrix(type = "Diag", nrow = 5, ncol = 5, free = TRUE , values = 1 , name = "U"), 
##' 	mxAlgebra(A %*% L %*% t(A) + U, name = "R"),
##' 	mxExpectationNormal(covariance = "R", dimnames = names(demoOneFactor)), 
##' 	mxFitFunctionML(), 
##' 	mxData(cov(demoOneFactor), type = "cov", numObs = 500), 
##' 	mxComputeSequence(
##' 		list(mxComputeNumericDeriv(), mxComputeReportDeriv())
##' 	)
##' )
##' factorModelFit <- mxRun(factorModel)
##' factorModelFit$output$hessian

mxComputeNumericDeriv <- function(freeSet=NA_character_, ..., fitfunction='fitfunction',
				      parallel=TRUE, stepSize=0.0001, iterations=4L, verbose=0L)
{
	garbageArguments <- list(...)
	if (length(garbageArguments) > 0) {
		stop("mxComputeNumericDeriv does not accept values for the '...' argument")
	}

	new("MxComputeNumericDeriv", freeSet, fitfunction, parallel, stepSize, iterations, verbose)
}

setMethod("displayCompute", signature("MxComputeNumericDeriv"),
	  function(Ob, indent) {
		  callNextMethod();
		  sp <- paste(rep('  ', indent), collapse="")
		  for (sl in c("fitfunction", "parallel", "stepSize", "iterations", "verbose")) {
			  slname <- paste("$", sl, sep="")
			  if (is.character(slot(Ob, sl))) {
				  cat(sp, slname, ":", omxQuotes(slot(Ob, sl)), '\n')
			  } else {
				  cat(sp, slname, ":", slot(Ob, sl), '\n')
			  }
		  }
		  invisible(Ob)
	  })

#----------------------------------------------------

setClass(Class = "MxComputeStandardError",
	 contains = "BaseCompute")

setMethod("initialize", "MxComputeStandardError",
	  function(.Object, freeSet) {
		  .Object@name <- 'compute'
		  .Object@freeSet <- freeSet
		  .Object
	  })

##' Compute standard errors given the Hessian or inverse Hessian
##'
##' @param freeSet names of matrices containing free variables
##' @aliases
##' MxComputeStandardError-class

mxComputeStandardError <- function(freeSet=NA_character_) {
	new("MxComputeStandardError", freeSet)
}

#----------------------------------------------------

setClass(Class = "MxComputeHessianQuality",
	 contains = "BaseCompute")

setMethod("initialize", "MxComputeHessianQuality",
	  function(.Object, freeSet) {
		  .Object@name <- 'compute'
		  .Object@freeSet <- freeSet
		  .Object
	  })

##' Compute the quality of the Hessian
##'
##' Tests whether the Hessian is positive definite
##' (model$output$infoDefinite) and, if so, computes the condition
##' number (model$output$conditionNumber). See Luenberger & Ye (2008)
##' Second Order Test (p. 190) and Condition Number (p. 239).
##' 
##' @param freeSet names of matrices containing free variables
##' @aliases
##' MxComputeHessianQuality-class
##' @references
##' Luenberger, D. G. & Ye, Y. (2008). Linear and nonlinear programming. Springer.

mxComputeHessianQuality <- function(freeSet=NA_character_) {
	new("MxComputeHessianQuality", freeSet)
}

#----------------------------------------------------

setClass(Class = "MxComputeReportDeriv",
	 contains = "BaseCompute")

setMethod("initialize", "MxComputeReportDeriv",
	  function(.Object, freeSet) {
		  .Object@name <- 'compute'
		  .Object@freeSet <- freeSet
		  .Object
	  })

##' Report derivatives
##'
##' Copy the internal gradient and Hessian back to R.
##'
##' @param freeSet names of matrices containing free variables
##' @aliases
##' MxComputeReportDeriv-class

mxComputeReportDeriv <- function(freeSet=NA_character_) {
	new("MxComputeReportDeriv", freeSet)
}

#----------------------------------------------------

setClass(Class = "MxComputeSequence",
	 contains = "ComputeSteps")

setMethod("initialize", "MxComputeSequence",
	  function(.Object, steps, freeSet) {
		  .Object@name <- 'compute'
		  .Object@steps <- steps
		  .Object@freeSet <- freeSet
		  .Object
	  })

##' Invoke a series of compute objects in sequence
##'
##' @param steps a list of compute objects
##' @param ... Not used; forces argument 'freeSet' to be specified by name.
##' @param freeSet Names of matrices containing free parameters.
##' @aliases
##' MxComputeSequence-class
mxComputeSequence <- function(steps=list(), ..., freeSet=NA_character_) {
	garbageArguments <- list(...)
	if (length(garbageArguments) > 0) {
		stop("mxComputeSequence does not accept values for the '...' argument")
	}

	new("MxComputeSequence", steps=steps, freeSet)
}

##' Compute nothing
##'
##' Note that this compute plan actually does nothing whereas
##' \code{mxComputeOnce("expectation", "nothing")} may remove the
##' prediction of an expectation.
##' 
mxComputeNothing <- function() {
	mxComputeSequence(freeSet=c())
}

setMethod("displayCompute", signature("MxComputeSequence"),
	  function(Ob, indent) {
		  callNextMethod();
		  sp <- paste(rep('  ', indent), collapse="")
		  for (step in 1:length(Ob@steps)) {
			  cat(sp, "steps[[", step, "]] :", '\n')
			  displayCompute(Ob@steps[[step]], indent+1L)
		  }
		  invisible(Ob)
	  })

convertComputes <- function(flatModel, model) {
	if (is.null(flatModel@compute)) return()
	convertForBackend(flatModel@compute, flatModel, model)
}

updateModelCompute <- function(model, computes) {
	if (is.null(model@compute)) return()
	updateFromBackend(model@compute, computes)
}

##' Sparse symmetric matrix invert
##'
##' This API is visible to permit testing. Please do not use.
##'
##' @param mat the matrix to invert
imxSparseInvert <- function(mat) .Call(sparseInvert_wrapper, mat)
