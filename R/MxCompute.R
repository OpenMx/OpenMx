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

##' BaseCompute
##'
##' This is an internal class and should not be used directly.
##'
##' @aliases
##' $,BaseCompute-method
##' $<-,BaseCompute-method
##' print,BaseCompute-method
##' show,BaseCompute-method
##' @seealso
##' \link{mxComputeEM}, \link{mxComputeGradientDescent}, \link{mxComputeHessianQuality},
##' \link{mxComputeIterate}, \link{mxComputeNewtonRaphson}, \link{mxComputeNumericDeriv}
##' @rdname BaseCompute-class
setClass(Class = "BaseCompute",
	 representation = representation(
	   id = "integer",
	     freeSet = "MxOptionalChar",
	     output = "list",
	     debug = "list",
	     .persist = "logical",
	   "VIRTUAL"),
	 contains = "MxBaseNamed")

##' @title MxCompute
##' @name MxCompute-class
##'
##' @description
##' This is an internal class and should not be used directly.
##'
##' @aliases
##' MxCompute
##' MxCompute-class
##' @rdname MxCompute-class
setClassUnion("MxCompute", c("NULL", "BaseCompute"))

setGeneric("displayCompute",
	   function(Ob, indent) {
		   return(standardGeneric("displayCompute"))
	   })

setMethod("displayCompute", signature(Ob="BaseCompute", indent="integer"),
	  function(Ob, indent) {
		  sp <- paste(rep('  ', indent), collapse="")
		  cat(sp, class(Ob), omxQuotes(Ob@name), '\n')
#		  cat(sp, "$id :", Ob@id, '\n')   # only of interest to developers and introduces visual noise
		  cat(sp, "$freeSet :", omxQuotes(Ob@freeSet), '\n')
		  if (length(Ob$output)) {
			  for (elem in names(Ob$output)) {
				  stuff <- Ob@output[[elem]]
				  if (is.list(stuff)) {
					  cat(sp, "$output[[", omxQuotes(elem), "]] : ...", '\n')
				  } else {
					  cat(sp, "$output[[", omxQuotes(elem), "]] :", stuff, '\n')
				  }
			  }
		  }
		  if (length(Ob$debug)) {
			  for (elem in names(Ob$debug)) {
				  cat(sp, "$debug[[", omxQuotes(elem), "]] : ...", '\n')
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

setMethod("names", "BaseCompute", slotNames)

#----------------------------------------------------

setClass(Class = "MxComputeOnce",
	 contains = "BaseCompute",
	 representation = representation(
	     from = "MxCharOrNumber",
	     what = "MxOptionalChar",
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
				stop(paste("Can only apply MxComputeOnce to MxFitFunction or MxExpectation not",
					   deparse(.Object@from)))
			}
			if (!any(is.na(expNum))) {
					# Usually negative numbers indicate matrices; not here
				.Object@from <- - expNum
			} else {
				.Object@from <- algNum - 1L
			}
		}
		.Object
	})

setMethod("initialize", "MxComputeOnce",
	  function(.Object, from, what, how, freeSet, verbose, .is.bestfit) {
		  .Object@name <- 'compute'
		  .Object@.persist <- TRUE
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
##' @param from the object to perform the computation (a vector of expectation or fit function names)
##' @param what what to compute
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

mxComputeOnce <- function(from, what=NULL, how=NULL, ...,
			  freeSet=NA_character_, verbose=0L, .is.bestfit=FALSE) {
	garbageArguments <- list(...)
	if (length(garbageArguments) > 0) {
		stop("mxComputeOnce does not accept values for the '...' argument")
	}
	if (length(from) == 0) warning("mxComputeOnce from nothing will have no effect")
	verbose <- as.integer(verbose)
	new("MxComputeOnce", from, what, how, freeSet, verbose, .is.bestfit)
}

setMethod("displayCompute", signature(Ob="MxComputeOnce", indent="integer"),
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
	     availableEngines = "character",
	     tolerance = "numeric",
	   nudgeZeroStarts = "MxCharOrLogical",
	   .excludeVars = "MxOptionalChar",
	   verbose = "integer",
	     maxMajorIter = "integer",
	     gradientAlgo = "character",
	     gradientIterations = "integer",
	   gradientStepSize = "numeric",
	   defaultCImethod = "character",
	     warmStart = "MxOptionalMatrix"))  # rename to 'preconditioner'?

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
	  function(.Object, freeSet, engine, fit, useGradient, verbose, tolerance, warmStart,
		   nudgeZeroStarts, maxMajorIter, gradientAlgo, gradientIterations, gradientStepSize) {
		  .Object@name <- 'compute'
		  .Object@.persist <- TRUE
		  .Object@freeSet <- freeSet
		  .Object@fitfunction <- fit
		  .Object@engine <- engine
		  .Object@defaultCImethod <- 'none'
		  if (engine == 'SLSQP') .Object@defaultCImethod <- 'ineq'
		  .Object@useGradient <- useGradient
		  .Object@verbose <- verbose
		  .Object@tolerance <- tolerance
		  .Object@warmStart <- warmStart
		  .Object@nudgeZeroStarts <- nudgeZeroStarts
		  .Object@maxMajorIter <- maxMajorIter
		  .Object@gradientAlgo <- gradientAlgo
		  .Object@gradientIterations <- gradientIterations
		  .Object@gradientStepSize <- gradientStepSize
		  .Object@availableEngines <- c("CSOLNP", "SLSQP")
		  if (imxHasNPSOL()) {
			  .Object@availableEngines <- c(.Object@availableEngines, "NPSOL")
		  }
		  .Object@.excludeVars <- c()
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
##' function. The fully open-source CRAN version of OpenMx offers 2 choices,
##' SLSQP (from the NLOPT collection) and CSOLNP.  The OpenMx Team's version of
##' OpenMx offers the choice of three optimizers: SLSQP, CSOLNP, and NPSOL.
##'
##' One of the most important options for SLSQP is
##' \code{gradientAlgo}. By default, the \code{central} method 
##' is used.  This method requires 2 times 
##' \code{gradientIterations} function evaluations per parameter 
##' per gradient.  The \code{central} method can be much more accurate 
##' than the \code{forward} method, which requires 
##' \code{gradientIterations} function evaluations per parameter per 
##' gradient.  The \code{forward} method is faster, and often works
##' well enough, but can result in imprecise gradient estimations that
##' may not allow SLSQP to fully optimize a given model, possibly 
##' resulting in code RED (status code 5 or 6).
##' 
##' Currently, only SLSQP uses arguments \code{gradientIterations} 
##' and \code{gradientAlgo}.  CSOLNP always uses the \code{forward} 
##' method; NPSOL usually uses the \code{forward} method, but 
##' adaptively switches to \code{central} under certain circumstances.
##' 
##' CSOLNP uses the value of argument \code{gradientStepSize} as-is, 
##' whereas SLSQP internally scales it by a factor of 100. The
##' purpose of this transformation is to obtain roughly the same
##' accuracy given other differences in numerical procedure.
##' NPSOL ignores \code{gradientStepSize}, and instead uses a function
##' of \link{mxOption} \dQuote{Function precision} to determine its gradient
##' step size.
##' 
##' Currently, only SLSQP and NPSOL can use analytic gradients, 
##' and only NPSOL uses \code{warmStart}.
##'
##' @param freeSet names of matrices containing free parameters.
##' @param ...  Not used.  Forces remaining arguments to be specified by name.
##' @param engine specific 'NPSOL', 'SLSQP', or 'CSOLNP'
##' @param fitfunction name of the fitfunction (defaults to 'fitfunction')
##' @param verbose level of debugging output
##' @param tolerance how close to the optimum is close enough (also known as the optimality tolerance)
##' @param useGradient whether to use the analytic gradient (if available)
##' @param warmStart a Cholesky factored Hessian to use as the NPSOL Hessian starting value (preconditioner)
##' @param nudgeZeroStarts whether to nudge any zero starting values prior to optimization (default TRUE)
##' @param maxMajorIter maximum number of major iterations
##' @param gradientAlgo one of c('forward','central')
##' @param gradientIterations number of Richardson iterations to use for the gradient
##' @param gradientStepSize the step size for the gradient
##' @aliases
##' MxComputeGradientDescent-class
##' @references
##' Luenberger, D. G. & Ye, Y. (2008). \emph{Linear and nonlinear programming.} Springer.
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
				     tolerance=NA_real_, useGradient=NULL, warmStart=NULL,
				     nudgeZeroStarts=mxOption(NULL,"Nudge zero starts"), maxMajorIter=NULL,
				     gradientAlgo=mxOption(NULL, "Gradient algorithm"),
				     gradientIterations=imxAutoOptionValue("Gradient iterations"),
				     gradientStepSize=imxAutoOptionValue("Gradient step size")) {

	garbageArguments <- list(...)
	if (length(garbageArguments) > 0) {
		stop("mxComputeGradientDescent does not accept values for the '...' argument")
	}
	if (missing(engine)) {
		engine <- options()$mxOptions[["Default optimizer"]]
	}

	if (!is.null(warmStart) && engine != "NPSOL") {
		stop("Only NPSOL supports warmStart")
	}
	verbose <- as.integer(verbose)
	maxMajorIter <- as.integer(maxMajorIter)
	gradientIterations <- as.integer(gradientIterations)

	new("MxComputeGradientDescent", freeSet, engine, fitfunction, useGradient, verbose,
	    tolerance, warmStart, nudgeZeroStarts, maxMajorIter,
	    gradientAlgo, gradientIterations, gradientStepSize)
}

setMethod("displayCompute", signature(Ob="MxComputeGradientDescent", indent="integer"),
	  function(Ob, indent) {
		  callNextMethod();
		  sp <- paste(rep('  ', indent), collapse="")
		  for (sl in c("engine", "fitfunction", "verbose", "tolerance", "useGradient",
			       "nudgeZeroStarts", "maxMajorIter",
			       "gradientAlgo", "gradientIterations", "gradientStepSize")) {
			  val <- slot(Ob, sl)
			  if (length(val)==0 || is.na(val)) next
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

setClass(Class = "MxComputeTryHard",
	 contains = "BaseCompute",
	 representation = representation(
	     plan = "MxCompute",
	     verbose = "integer",
	     location = "numeric",
	     scale = "numeric",
	     maxRetries = "integer"))

setMethod("assignId", signature("MxComputeTryHard"),
	function(.Object, id, defaultFreeSet) {
		.Object <- callNextMethod()
		defaultFreeSet <- .Object@freeSet
		id <- .Object@id
		for (sl in c('plan')) {
			slot(.Object, sl) <- assignId(slot(.Object, sl), id, defaultFreeSet)
			id <- slot(.Object, sl)@id + 1L
		}
		.Object@id <- id
		.Object
	})

setMethod("getFreeVarGroup", signature("MxComputeTryHard"),
	function(.Object) {
		result <- callNextMethod()
		for (step in c(.Object@plan)) {
			got <- getFreeVarGroup(step)
			if (length(got)) result <- append(result, got)
		}
		result
	})

setMethod("qualifyNames", signature("MxComputeTryHard"),
	function(.Object, modelname, namespace) {
		.Object <- callNextMethod()
		for (sl in c('plan')) {
			slot(.Object, sl) <- qualifyNames(slot(.Object, sl), modelname, namespace)
		}
		.Object
	})

setMethod("convertForBackend", signature("MxComputeTryHard"),
	function(.Object, flatModel, model) {
		name <- .Object@name
		for (sl in c('plan')) {
			slot(.Object, sl) <- convertForBackend(slot(.Object, sl), flatModel, model)
		}
		.Object
	})

setMethod("updateFromBackend", signature("MxComputeTryHard"),
	function(.Object, computes) {
		.Object <- callNextMethod()
		for (sl in c('plan')) {
			slot(.Object, sl) <- updateFromBackend(slot(.Object, sl), computes)
		}
		.Object
	})

setMethod("initialize", "MxComputeTryHard",
	  function(.Object, freeSet, plan, verbose, location, scale, maxRetries) {
		  .Object@name <- 'compute'
		  .Object@.persist <- TRUE
		  .Object@freeSet <- freeSet
		  .Object@plan <- plan
		  .Object@verbose <- verbose
		  .Object@location <- location
		  .Object@scale <- scale
		  .Object@maxRetries <- maxRetries
		  .Object
	  })

##' Repeatedly attempt a compute plan until successful
##'
##' The provided compute plan is run until the status code indicates
##' success (0 or 1). It gives up after a small number of retries.
##'
##' Upon failure, start values are randomly perturbed.  Currently only
##' the uniform distribution is implemented.  The distribution is
##' parametrized by arguments \code{location} and \code{scale}.  The
##' location parameter is the distribution's median.  For the uniform
##' distribution, \code{scale} is the absolute difference between its
##' median and extrema (i.e., half the width of the rectangle).  Each
##' start value is multiplied by a random draw and then added to a
##' random draw from a distribution with the same \code{scale} but
##' with a median of zero.
##' 
##' @param plan compute plan to optimize the model
##' @param ...  Not used.  Forces remaining arguments to be specified by name.
##' @param freeSet names of matrices containing free variables
##' @param verbose level of debugging output
##' @param location location of the perturbation distribution
##' @param scale scale of the perturbation distribution
##' @param maxRetries maximum number of plan evaluations per invocation (including the first evaluation)
##' @seealso
##' \code{\link{mxTryHard}}
##' @aliases
##' MxComputeTryHard-class
##' @references
##' Shanno, D. F. (1985). On Broyden-Fletcher-Goldfarb-Shanno method. \emph{Journal of
##' Optimization Theory and Applications, 46}(1), 87-94.
mxComputeTryHard <- function(plan, ..., freeSet=NA_character_, verbose=0L,
			     location=1.0, scale=0.25, maxRetries=3L)
{
	garbageArguments <- list(...)
	if (length(garbageArguments) > 0) {
		stop("mxComputeTryHard does not accept values for the '...' argument")
	}
	verbose <- as.integer(verbose)
	maxRetries <- as.integer(maxRetries)
	new("MxComputeTryHard", freeSet, plan, verbose, location, scale, maxRetries)
}

setMethod("displayCompute", signature(Ob="MxComputeTryHard", indent="integer"),
	  function(Ob, indent) {
		  callNextMethod()
		  sp <- paste(rep('  ', indent), collapse="")
		  cat(sp, "$plan :", '\n')
		  displayCompute(Ob@plan, indent+1L)
		  for (sl in c("verbose","location","scale",'maxRetries')) {
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

setClass(Class = "MxComputeConfidenceInterval",
	 contains = "BaseCompute",
	 representation = representation(
	     plan = "MxCompute",
	   fitfunction = "MxCharOrNumber",
	     constraintType = "character",
	     verbose = "integer"))

setMethod("assignId", signature("MxComputeConfidenceInterval"),
	function(.Object, id, defaultFreeSet) {
		.Object <- callNextMethod()
		defaultFreeSet <- .Object@freeSet
		id <- .Object@id
		for (sl in c('plan')) {
			slot(.Object, sl) <- assignId(slot(.Object, sl), id, defaultFreeSet)
			id <- slot(.Object, sl)@id + 1L
		}
		.Object@id <- id
		.Object
	})

setMethod("getFreeVarGroup", signature("MxComputeConfidenceInterval"),
	function(.Object) {
		result <- callNextMethod()
		for (step in c(.Object@plan)) {
			got <- getFreeVarGroup(step)
			if (length(got)) result <- append(result, got)
		}
		result
	})

setMethod("qualifyNames", signature("MxComputeConfidenceInterval"),
	function(.Object, modelname, namespace) {
		.Object <- callNextMethod()
		for (sl in c('plan')) {
			slot(.Object, sl) <- qualifyNames(slot(.Object, sl), modelname, namespace)
		}
		for (sl in c('fitfunction')) {
			slot(.Object, sl) <- imxConvertIdentifier(slot(.Object, sl), modelname, namespace)
		}
		.Object
	})

setMethod("convertForBackend", signature("MxComputeConfidenceInterval"),
	function(.Object, flatModel, model) {
		name <- .Object@name
		for (sl in c('plan')) {
			slot(.Object, sl) <- convertForBackend(slot(.Object, sl), flatModel, model)
		}
		if (is.character(.Object@fitfunction)) {
			.Object@fitfunction <- imxLocateIndex(flatModel, .Object@fitfunction, .Object)
		}
		.Object
	})

setMethod("initialize", "MxComputeConfidenceInterval",
	  function(.Object, freeSet, plan, verbose, fitfunction, constraintType) {
		  .Object@name <- 'compute'
		  .Object@.persist <- TRUE
		  .Object@freeSet <- freeSet
		  .Object@plan <- plan
		  .Object@verbose <- verbose
		  .Object@fitfunction <- fitfunction
		  .Object@constraintType <- constraintType
		  .Object
	  })

##' Find likelihood-based confidence intervals
##'
##' There are various equivalent ways to pose the optimization
##' problems required to estimate confindence intervals. Most accurate
##' solutions are achieved when the problem is posed using non-linear
##' constraints. However, the available optimizers (NPSOL, SLSQP, and
##' CSOLNP) often have difficulty with non-linear
##' constraints.
##' 
##' @param plan compute plan to optimize the model
##' @param ...  Not used.  Forces remaining arguments to be specified by name.
##' @param freeSet names of matrices containing free variables
##' @param verbose level of debugging output
##' @param engine deprecated
##' @param fitfunction the name of the deviance function
##' @param tolerance deprecated
##' @param constraintType one of c('ineq', 'none')
##' @references
##' Neale, M. C. & Miller M. B. (1997). The use of likelihood based
##' confidence intervals in genetic models.  \emph{Behavior Genetics,
##' 27}(2), 113-120.
##' 
##' Pek, J. & Wu, H. (2015). Profile likelihood-based confidence intervals and regions for structural equation models.
##' \emph{Psychometrica, 80}(4), 1123-1145.
##'
##' Wu, H. & Neale, M. C. (2012). Adjusted confidence intervals for a
##' bounded parameter. \emph{Behavior genetics, 42}(6), 886-898.
##' @aliases
##' MxComputeConfidenceInterval-class

mxComputeConfidenceInterval <- function(plan, ..., freeSet=NA_character_, verbose=0L,
					engine=NULL, fitfunction='fitfunction',
					tolerance=NA_real_, constraintType='none') {

	garbageArguments <- list(...)
	if (length(garbageArguments) > 0) {
		stop("mxComputeConfidenceInterval does not accept values for the '...' argument")
	}
	verbose <- as.integer(verbose)
	new("MxComputeConfidenceInterval", freeSet, plan, verbose, fitfunction, constraintType)
}

setMethod("updateFromBackend", signature("MxComputeConfidenceInterval"),
	function(.Object, computes) {
		.Object <- callNextMethod()
		for (sl in c('plan')) {
			slot(.Object, sl) <- updateFromBackend(slot(.Object, sl), computes)
		}
		.Object
	})

setMethod("displayCompute", signature(Ob="MxComputeConfidenceInterval", indent="integer"),
	  function(Ob, indent) {
		  callNextMethod()
		  sp <- paste(rep('  ', indent), collapse="")
		  cat(sp, "$plan :", '\n')
		  displayCompute(Ob@plan, indent+1L)
		  for (sl in c("verbose")) {
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
		  .Object@.persist <- TRUE
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
##' fit function. Box constraints are supported. Parameters can
##' approach box constraints but will not leave the feasible region
##' (even by some small epsilon>0). Non-finite fit values are
##' interpreted as soft feasibility constraints. That is, when a
##' non-finite fit is encountered, line search is continued after the
##' step size is multipled by 10%. Comprehensive diagnostics are
##' available by increasing the verbose level.
##'
##' @param freeSet names of matrices containing free variables
##' @param ...  Not used.  Forces remaining arguments to be specified by name.
##' @param fitfunction name of the fitfunction (defaults to 'fitfunction')
##' @param maxIter maximum number of iterations
##' @param tolerance optimization is considered converged when the maximum relative change in fit is less than tolerance
##' @param verbose level of debugging output
##' @aliases
##' MxComputeNewtonRaphson-class
##' @references
##' Luenberger, D. G. & Ye, Y. (2008). \emph{Linear and nonlinear programming.} Springer.

mxComputeNewtonRaphson <- function(freeSet=NA_character_, ..., fitfunction='fitfunction', maxIter = 100L,
				   tolerance=1e-12, verbose=0L)
{
	garbageArguments <- list(...)
	if (length(garbageArguments) > 0) {
		stop("mxComputeNewtonRaphson does not accept values for the '...' argument")
	}

	verbose <- as.integer(verbose)
	maxIter <- as.integer(maxIter)
	new("MxComputeNewtonRaphson", freeSet, fitfunction, maxIter, tolerance, verbose)
}

setMethod("displayCompute", signature(Ob="MxComputeNewtonRaphson", indent="integer"),
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
	   verbose = "integer",
	   maxDuration = "numeric"))

setMethod("initialize", "MxComputeIterate",
	  function(.Object, steps, maxIter, tolerance, verbose, freeSet, maxDuration) {
		  .Object@name <- 'compute'
		  .Object@.persist <- TRUE
		  .Object@steps <- steps
		  .Object@maxIter <- maxIter
		  .Object@tolerance <- tolerance
		  .Object@verbose <- verbose
		  .Object@freeSet <- freeSet
		  .Object@maxDuration <- maxDuration
		  .Object
	  })

##' Repeatedly invoke a series of compute objects until change is less than tolerance
##'
##' One step (typically the last) must compute the fit or maxAbsChange.
##'
##' @param steps a list of compute objects
##' @param ...  Not used.  Forces remaining arguments to be specified by name.
##' @param maxIter the maximum number of iterations
##' @param tolerance iterates until maximum relative change is less than tolerance
##' @param verbose level of debugging output
##' @param freeSet Names of matrices containing free variables.
##' @param maxDuration the maximum amount of time (in seconds) to iterate
##' @aliases
##' MxComputeIterate-class
mxComputeIterate <- function(steps, ..., maxIter=500L, tolerance=1e-9, verbose=0L, freeSet=NA_character_,
			     maxDuration=as.numeric(NA)) {
	garbageArguments <- list(...)
	if (length(garbageArguments) > 0) {
		stop("mxComputeIterate does not accept values for the '...' argument")
	}

	verbose <- as.integer(verbose)
	maxIter <- as.integer(maxIter)
	new("MxComputeIterate", steps=steps, maxIter=maxIter, tolerance=tolerance,
	    verbose, freeSet, maxDuration)
}

setMethod("displayCompute", signature(Ob="MxComputeIterate", indent="integer"),
	  function(Ob, indent) {
		  callNextMethod();
		  sp <- paste(rep('  ', indent), collapse="")
		  cat(sp, "maxIter :", Ob@maxIter, '\n')
		  cat(sp, "maxDuration :", Ob@maxDuration, '\n')
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
	     estep = "MxCompute",
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
		.Object <- callNextMethod()
		defaultFreeSet <- .Object@freeSet
		id <- .Object@id
		for (sl in c('estep', 'mstep')) {
			slot(.Object, sl) <- assignId(slot(.Object, sl), id, defaultFreeSet)
			id <- slot(.Object, sl)@id + 1L
		}
		.Object@id <- id 
		.Object
	})

setMethod("getFreeVarGroup", signature("MxComputeEM"),
	function(.Object) {
		result <- callNextMethod()
		for (step in c(.Object@estep, .Object@mstep)) {
			got <- getFreeVarGroup(step)
			if (length(got)) result <- append(result, got)
		}
		result
	})

setMethod("qualifyNames", signature("MxComputeEM"),
	function(.Object, modelname, namespace) {
		.Object <- callNextMethod()
		for (sl in c('observedFit')) {
			slot(.Object, sl) <- imxConvertIdentifier(slot(.Object, sl), modelname, namespace)
		}
		for (sl in c('estep', 'mstep')) {
			slot(.Object, sl) <- qualifyNames(slot(.Object, sl), modelname, namespace)
		}
		.Object@infoArgs$fitfunction <-
		    imxConvertIdentifier(.Object@infoArgs$fitfunction, modelname, namespace)
		.Object
	})

setMethod("convertForBackend", signature("MxComputeEM"),
	function(.Object, flatModel, model) {
		name <- .Object@name
		if (length(.Object@observedFit) != 1) stop("MxComputeEM requires a single observedFit function")
		if (any(!is.integer(.Object@observedFit))) {
			algNum <- match(.Object@observedFit, append(names(flatModel@algebras),
								    names(flatModel@fitfunctions)))
			if (any(is.na(algNum))) {
				stop(paste("MxComputeEM: observedFit fit function", .Object@observedFit, "not found"))
			}
			.Object@observedFit <- algNum - 1L
		}
		for (sl in c('estep', 'mstep')) {
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
		for (sl in c('estep', 'mstep')) {
			slot(.Object, sl) <- updateFromBackend(slot(.Object, sl), computes)
		}
		.Object
	})

setMethod("initialize", "MxComputeEM",
	  function(.Object, estep, mstep, observedFit, maxIter, tolerance,
		   verbose, accel, information, freeSet, infoArgs) {
		  .Object@name <- 'compute'
		  .Object@.persist <- TRUE
		  .Object@estep <- estep
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
##' The EM algorithm does not produce a parameter covariance matrix
##' for standard errors. The Oakes (1999) direct method and S-EM, an
##' implementation of Meng & Rubin (1991), are included.
##'
##' Ramsay (1975) was recommended in Bock, Gibbons, & Muraki (1988).
##'
##' @param expectation a vector of expectation names (DEPRECATED)
##' @param predict what to predict from the observed data (DEPRECATED)
##' @param mstep a compute plan to optimize the completed data model
##' @param observedFit the name of the observed data fit function (defaults to "fitfunction")
##' @param ...  Not used.  Forces remaining arguments to be specified by name.
##' @param maxIter maximum number of iterations
##' @param tolerance optimization is considered converged when the maximum relative change in fit is less than tolerance
##' @param verbose level of diagnostic output
##' @param freeSet names of matrices containing free variables
##' @param accel name of acceleration method ("varadhan2008" or "ramsay1975")
##' @param information name of information matrix approximation method
##' @param infoArgs arguments to control the information matrix method
##' @param estep a compute plan to perform the expectation step
##' @aliases
##' MxComputeEM-class
##' @references
##'
##' Bock, R. D., Gibbons, R., & Muraki, E. (1988). Full-information
##' item factor analysis. \emph{Applied Psychological Measurement,
##' 6}(4), 431-444.
##' 
##' Dempster, A. P., Laird, N. M., & Rubin, D. B. (1977). Maximum likelihood from
##' incomplete data via the EM algorithm. \emph{Journal of the Royal Statistical Society.
##' Series B (Methodological)}, 1-38.
##'
##' Meng, X.-L. & Rubin, D. B. (1991). Using EM to obtain asymptotic variance-covariance
##' matrices: The SEM algorithm. \emph{Journal of the American Statistical Association,
##' 86} (416), 899-909.
##'
##' Oakes, D. (1999). Direct calculation of the information matrix via
##' the EM algorithm.  \emph{Journal of the Royal Statistical Society:
##' Series B (Statistical Methodology), 61}(2), 479-482.
##' 
##' Ramsay, J. O. (1975). Solving implicit equations in psychometric data analysis.
##' \emph{Psychometrika, 40} (3), 337-360.
##'
##' Varadhan, R. & Roland, C. (2008). Simple and globally convergent
##' methods for accelerating the convergence of any EM
##' algorithm. \emph{Scandinavian Journal of Statistics, 35}, 335-353.
mxComputeEM <- function(expectation=NULL, predict=NA_character_, mstep, observedFit="fitfunction", ...,
			maxIter=500L, tolerance=1e-9, verbose=0L, freeSet=NA_character_,
			accel="varadhan2008", information=NA_character_, infoArgs=list(), estep=NULL) {
	garbageArguments <- list(...)
	if (length(garbageArguments) > 0) {
		stop("mxComputeEM does not accept values for the '...' argument")
	}
	verbose <- as.integer(verbose)
	maxIter <- as.integer(maxIter)
	accel <- as.character(accel)

	# backward compatibility for original API
	if (length(expectation)) {
		if (length(estep)) stop("You cannot provide both 'expectation' and 'estep' arguments")
		estep <- mxComputeOnce(expectation, predict)
		mstep <- mxComputeSequence(list(mstep, mxComputeOnce(expectation)))
	}

	new("MxComputeEM", estep, mstep, observedFit, maxIter=maxIter,
	    tolerance=tolerance, verbose, accel, information, freeSet, infoArgs)
}

setMethod("displayCompute", signature(Ob="MxComputeEM", indent="integer"),
	  function(Ob, indent) {
		  callNextMethod();
		  sp <- paste(rep('  ', indent), collapse="")
		  cat(sp, "$estep :", '\n')
		  displayCompute(Ob@estep, indent+1L)
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
#Mike Hunter's "better match.arg" function
match.barg <- function (arg, choices, several.ok = FALSE) 
{
	if (missing(choices)) {
		formal.args <- formals(sys.function(sys.parent()))
		choices <- eval(formal.args[[deparse(substitute(arg))]])
	}
	if (is.null(arg)) 
		return(choices[1L])
	else if (!is.character(arg)) 
		stop("'arg' must be NULL or a character vector")
	if (!several.ok) {
		if (identical(arg, choices)) 
			return(arg[1L])
		if (length(arg) > 1L) 
			stop("'arg' must be of length 1")
	}
	else if (length(arg) == 0L) 
		stop("'arg' must be of length >= 1")
	i <- pmatch(arg, choices, nomatch = 0L, duplicates.ok = TRUE)
	if (all(i == 0L)) 
		stop(gettextf("%s should be one of %s", omxQuotes(arg), omxQuotes(choices)
		), domain = NA)
	i <- i[i > 0L]
	if (!several.ok && length(i) > 1) 
		stop("there is more than one match in 'match.arg'")
	choices[i]
}

mxComputeNelderMead <- function(
	freeSet=NA_character_, fitfunction="fitfunction", verbose=0L, 
	nudgeZeroStarts=mxOption(NULL,"Nudge zero starts"), 
	maxIter=NULL,	...,
	alpha=1, betao=0.5, betai=0.5, gamma=2, sigma=0.5, bignum=1e35, 
	iniSimplexType=c("regular","right","smartRight","random"),
	iniSimplexEdge=1, iniSimplexMat=NULL, greedyMinimize=FALSE, 
	altContraction=FALSE, degenLimit=0, stagnCtrl=c(-1L,-1L),
	validationRestart=TRUE,
	xTolProx=1e-8, fTolProx=1e-8,
	doPseudoHessian=FALSE,
	ineqConstraintMthd=c("soft","eqMthd"), 
	eqConstraintMthd=c("GDsearch","soft","backtrack","l1p"),
	backtrackCtrl=c(0.5,5)
	){
	garbageArguments <- list(...)
	if (length(garbageArguments) > 0) {
		stop("mxComputeNelderMead() does not accept values for the '...' argument")
	}
	verbose <- as.integer(verbose[1])
	maxIter <- as.integer(maxIter[1])
	if(is.character(nudgeZeroStarts[1])){
		if(substr(nudgeZeroStarts[1],1,1) %in% c("Y","y")){nudgeZeroStarts <- TRUE}
		else if(substr(nudgeZeroStarts[1],1,1) %in% c("N","n")){nudgeZeroStarts <- FALSE}
		else{stop("unrecognized character string provided as argument 'nudgeZeroStarts'")}
	}
	alpha <- as.numeric(alpha[1])
	if(alpha<=0){stop("reflection coefficient 'alpha' must be positive")}
	betao <- as.numeric(betao[1])
	betai <- as.numeric(betai[1])
	if(any(betao<=0, betao>=1, betai<=0, betai>=1)){
		stop("contraction coefficients 'betao' and 'betai' must both be within unit interval (0,1)")
	}
	gamma <- as.numeric(gamma[1])
	#Allow user to provide non-positive gamma to "turn off" expanstion transformations:
	if(gamma>0 && gamma<=alpha){
		stop("if positive, expansion coefficient 'gamma' must be greater than reflection coefficient 'alpha'")
	}
	sigma <- as.numeric(sigma[1])
	#Allow user to provide non-positive sigma to "turn off" shrinks:
	if(sigma>=1){stop("shrink coefficient 'sigma' must be less than 1.0")}
	bignum <- as.numeric(bignum[1])
	iniSimplexType <- as.character(match.barg(iniSimplexType,c("regular","right","smartRight","random")))
	iniSimplexEdge <- as.numeric(iniSimplexEdge[1])
	if(length(iniSimplexMat)){
		iniSimplexMat <- as.matrix(iniSimplexMat)
		iniSimplexColnames <- colnames(iniSimplexMat)
	}
	else{
		iniSimplexColnames <- NULL
		iniSimplexMat <- NULL
	}
	greedyMinimize <- as.logical(greedyMinimize[1])
	altContraction <- as.logical(altContraction[1])
	degenLimit <- as.numeric(degenLimit[1])
	if(degenLimit<0 || degenLimit>pi){
		stop("'degenLimit' must be within interval [0,pi]")
	}
	if(length(stagnCtrl)<2){stop("'stagnCtrl' must be an integer vector of length 2")}
	stagnCtrl <- as.integer(stagnCtrl[1:2])
	validationRestart <- as.logical(validationRestart[1])
	xTolProx <- as.numeric(xTolProx[1])
	fTolProx <- as.numeric(fTolProx[1])
	backtrackCtrl=c(0.5,5)
	doPseudoHessian <- as.logical(doPseudoHessian[1])
	ineqConstraintMthd <- as.character(match.barg(ineqConstraintMthd,c("soft","eqMthd")))
	eqConstraintMthd <- as.character(match.barg(eqConstraintMthd,c("GDsearch","soft","backtrack","l1p")))
	if(length(backtrackCtrl)<2){stop("'backtrackCtrl' must be a numeric vector of length 2")}
	backtrackCtrl1 <- as.numeric(backtrackCtrl[1])
	backtrackCtrl2 <- as.integer(backtrackCtrl[2])
	return(new("MxComputeNelderMead", freeSet, fitfunction, verbose, nudgeZeroStarts, maxIter, alpha, 
						 betao, betai, gamma, sigma, bignum, iniSimplexType, iniSimplexEdge, iniSimplexMat, 
						 iniSimplexColnames, validationRestart,
						 greedyMinimize, altContraction, degenLimit, stagnCtrl, xTolProx, fTolProx, 
						 ineqConstraintMthd, eqConstraintMthd, backtrackCtrl1, backtrackCtrl2, doPseudoHessian))
}

setClass(
	Class="MxComputeNelderMead",
	contains="BaseCompute",
	representation=representation(
		fitfunction="MxCharOrNumber",
		verbose="integer",
		nudgeZeroStarts="MxCharOrLogical",
		maxIter="integer",
		defaultMaxIter="logical",
		.excludeVars="MxOptionalChar",
		alpha="numeric",
		betao="numeric",
		betai="numeric",
		gamma="numeric",
		sigma="numeric",
		bignum="numeric",
		iniSimplexType="character",
		iniSimplexEdge="numeric",
		iniSimplexMat="MxOptionalMatrix",
		.iniSimplexColnames="MxOptionalChar",
		greedyMinimize="logical",
		altContraction="logical",
		degenLimit="numeric",
		stagnCtrl="integer",
		validationRestart="logical",
		xTolProx="numeric",
		fTolProx="numeric",
		doPseudoHessian="logical",
		ineqConstraintMthd="character",
		eqConstraintMthd="character",
		backtrackCtrl1="numeric",
		backtrackCtrl2="integer"))

#TODO: a user or developer might someday want to directly use this low-level 'initialize' method instead of the high-level constructor function,
#so typecasting should also occur here:
setMethod(
	"initialize", "MxComputeNelderMead",
	function(.Object, freeSet, fitfunction, verbose, nudgeZeroStarts, maxIter, alpha, 
					 betao, betai, gamma, sigma, bignum, iniSimplexType, iniSimplexEdge, iniSimplexMat, 
					 iniSimplexColnames, validationRestart,
					 greedyMinimize, altContraction, degenLimit, stagnCtrl, xTolProx, fTolProx, 
					 ineqConstraintMthd, eqConstraintMthd, backtrackCtrl1, backtrackCtrl2, doPseudoHessian){
		.Object@name <- 'compute'
		.Object@.persist <- TRUE
		.Object@freeSet <- freeSet
		.Object@fitfunction <- fitfunction
		.Object@verbose <- verbose
		.Object@nudgeZeroStarts <- nudgeZeroStarts
		.Object@defaultMaxIter <- ifelse(length(maxIter),FALSE,TRUE)
		.Object@maxIter <- as.integer(maxIter)
		.Object@.excludeVars <- c()
		
		.Object@alpha <- alpha
		.Object@betao <- betao
		.Object@betai <- betai
		.Object@gamma <- gamma
		.Object@sigma <- sigma
		.Object@bignum <- bignum
		.Object@iniSimplexType <- iniSimplexType
		.Object@iniSimplexEdge <- iniSimplexEdge
		if(!length(iniSimplexMat)){.Object@iniSimplexMat <- NULL}
		else{.Object@iniSimplexMat <- iniSimplexMat}
		.Object@.iniSimplexColnames <- iniSimplexColnames
		.Object@greedyMinimize <- greedyMinimize
		.Object@altContraction <- altContraction
		.Object@degenLimit <- degenLimit
		.Object@stagnCtrl <- stagnCtrl
		.Object@validationRestart <- validationRestart
		.Object@xTolProx <- xTolProx
		.Object@fTolProx <- fTolProx
		.Object@doPseudoHessian <- doPseudoHessian
		.Object@ineqConstraintMthd <- ineqConstraintMthd
		.Object@eqConstraintMthd <- eqConstraintMthd
		.Object@backtrackCtrl1 <- backtrackCtrl1
		.Object@backtrackCtrl2 <- backtrackCtrl2
		.Object
	})

setMethod("qualifyNames", signature("MxComputeNelderMead"),
					function(.Object, modelname, namespace) {
						.Object <- callNextMethod()
						for (sl in c('fitfunction')) {
							slot(.Object, sl) <- imxConvertIdentifier(slot(.Object, sl), modelname, namespace)
						}
						if(length(.Object@iniSimplexMat)){.Object@.iniSimplexColnames <- colnames(.Object@iniSimplexMat)}
						.Object
					})

setMethod("convertForBackend", signature("MxComputeNelderMead"),
					function(.Object, flatModel, model) {
						name <- .Object@name
						if (is.character(.Object@fitfunction)) {
							.Object@fitfunction <- imxLocateIndex(flatModel, .Object@fitfunction, .Object)
						}
						.Object
					})

#----------------------------------------------------

setClass(Class = "MxComputeNumericDeriv",
	 contains = "BaseCompute",
	 representation = representation(
	   fitfunction = "MxCharOrNumber",
	     parallel = "logical",
	     stepSize = "numeric",
	     iterations = "integer",
	     verbose="integer",
	     knownHessian="MxOptionalMatrix",
	     checkGradient="logical",
	 hessian="logical"))

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
	  function(.Object, freeSet, fit, parallel, stepSize, iterations, verbose, knownHessian,
		   checkGradient, hessian) {
		  .Object@name <- 'compute'
		  .Object@.persist <- TRUE
		  .Object@freeSet <- freeSet
		  .Object@fitfunction <- fit
		  .Object@parallel <- parallel
		  .Object@stepSize <- stepSize
		  .Object@iterations <- iterations
		  .Object@verbose <- verbose
		  .Object@knownHessian <- knownHessian
		  .Object@checkGradient <- checkGradient
		  .Object@hessian <- hessian
		  .Object
	  })

adjustDefaultNumericDeriv <- function(m, iterations, stepSize) {
	for (nd in 1:length(m@compute@steps)) {
		if (is(m@compute@steps[[nd]], "MxComputeNumericDeriv")) {
			m@compute@steps[[nd]]$iterations <- iterations
			m@compute@steps[[nd]]$stepSize <- stepSize
			break
		}
	}
	m
}

##' Numerically estimate Hessian using Richardson extrapolation
##'
##' For N free parameters, Richardson extrapolation requires
##' (iterations * (N^2 + N)) function evaluations.
##' The implementation is closely based on the numDeriv R package.
##'
##' In addition to an estimate of the Hessian, forward, central, and
##' backward estimates of the gradient are made available in this
##' compute plan's output slot.
##' 
##' When \code{checkGradient=TRUE}, the central difference estimate of
##' the gradient is used to determine whether the first order
##' convergence criterion is met. In addition, the forward and
##' backward difference estimates of the gradient are compared for
##' symmetry. When sufficient asymmetry is detected, the standard
##' error is flagged. In the case, profile likelihood confidence
##' intervals should be used for inference instead of standard errors
##' (see \code{mxComputeConfidenceInterval}).
##'
##' @param freeSet names of matrices containing free variables
##' @param ...  Not used.  Forces remaining arguments to be specified by name.
##' @param fitfunction name of the fitfunction (defaults to 'fitfunction')
##' @param parallel whether to evaluate the fitfunction in parallel (defaults to TRUE)
##' @param stepSize starting set size (defaults to 0.0001)
##' @param iterations number of Richardson extrapolation iterations (defaults to 4L)
##' @param verbose Level of debugging output.
##' @param knownHessian an optional matrix of known Hessian entries
##' @param checkGradient whether to check the first order convergence criterion (gradient is near zero)
##' @param hessian whether to estimate the Hessian. If FALSE then only the gradient is estimated.
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
				  parallel=TRUE,
				  stepSize=imxAutoOptionValue("Gradient step size"),
				  iterations=4L, verbose=0L,
				  knownHessian=NULL, checkGradient=TRUE, hessian=TRUE)
{
	garbageArguments <- list(...)
	if (length(garbageArguments) > 0) {
		stop("mxComputeNumericDeriv does not accept values for the '...' argument")
	}

	verbose <- as.integer(verbose)
	iterations <- as.integer(iterations)

	if (!is.null(knownHessian)) {
		if (nrow(knownHessian) != ncol(knownHessian)) {
			stop("knownHessian must be square")
		}
		if (length(dimnames(knownHessian)) != 2 ||
		    any(rownames(knownHessian) != colnames(knownHessian))) {
			stop("knownHessian must have matching row and column names")
		}
		if (any(knownHessian != t(knownHessian))) {
			stop("knownHessian must be exactly symmetric")
		}
	}

	new("MxComputeNumericDeriv", freeSet, fitfunction, parallel, stepSize, iterations,
	    verbose, knownHessian, checkGradient, hessian)
}

setMethod("displayCompute", signature(Ob="MxComputeNumericDeriv", indent="integer"),
	  function(Ob, indent) {
		  callNextMethod();
		  sp <- paste(rep('  ', indent), collapse="")
		  for (sl in c("fitfunction", "parallel", "stepSize", "iterations",
			       "verbose", "knownHessian", 'checkGradient', 'hessian')) {
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
		  .Object@.persist <- TRUE
		  .Object@freeSet <- freeSet
		  .Object
	  })

##' Compute standard errors given the Hessian or inverse Hessian
##'
##' The fit is assumed to be in deviance units (-2 log likelihood).
##'
##' @param freeSet names of matrices containing free variables
##' @aliases
##' MxComputeStandardError-class

mxComputeStandardError <- function(freeSet=NA_character_) {
	new("MxComputeStandardError", freeSet)
}

#----------------------------------------------------

setClass(Class = "MxComputeHessianQuality",
	 contains = "BaseCompute",
	 representation = representation(
	     verbose = "integer"))

setMethod("initialize", "MxComputeHessianQuality",
	  function(.Object, freeSet, verbose) {
		  .Object@name <- 'compute'
		  .Object@.persist <- TRUE
		  .Object@freeSet <- freeSet
		  .Object@verbose <- verbose
		  .Object
	  })

##' Compute the quality of the Hessian
##'
##' Tests whether the Hessian is positive definite
##' (model$output$infoDefinite) and, if so, computes the approximate condition
##' number (model$output$conditionNumber). See Luenberger & Ye (2008)
##' Second Order Test (p. 190) and Condition Number (p. 239).
##'
##' The condition number is approximated by \eqn{\mathrm{norm}(H) *
##' \mathrm{norm}(H^{-1})}{norm(H) * norm(solve(H))} where H is the
##' Hessian. The norm is either the 1- or infinity-norm (both obtain
##' the same result due to symmetry).
##' 
##' @param freeSet names of matrices containing free variables
##' @param ...  Not used.  Forces remaining arguments to be specified by name.
##' @param verbose Level of debugging output.
##' @aliases
##' MxComputeHessianQuality-class
##' @references
##' Luenberger, D. G. & Ye, Y. (2008). Linear and nonlinear programming. Springer.

mxComputeHessianQuality <- function(freeSet=NA_character_, ..., verbose=0L) {
	garbageArguments <- list(...)
	if (length(garbageArguments) > 0) {
		stop("mxComputeHessianQuality does not accept values for the '...' argument")
	}

	new("MxComputeHessianQuality", freeSet, as.integer(verbose))
}

#----------------------------------------------------

setClass(Class = "MxComputeReportDeriv",
	 contains = "BaseCompute")

setMethod("initialize", "MxComputeReportDeriv",
	  function(.Object, freeSet) {
		  .Object@name <- 'compute'
		  .Object@.persist <- TRUE
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

setClass(Class = "MxComputeBootstrap",
	 contains = "BaseCompute",
	 representation = representation(
		 data = "MxCharOrNumber",
	     plan = "MxCompute",
	     replications = "integer",
	     verbose = "integer",
	     parallel = "logical",
	     OK = "character",
	     only = "integer"
	 ))

setMethod("initialize", "MxComputeBootstrap",
	  function(.Object, freeSet, data, plan, replications,
		   verbose, parallel, OK, only) {
		  .Object@name <- 'compute'
		  .Object@.persist <- TRUE
		  .Object@freeSet <- freeSet
		  .Object@data <- data
		  .Object@plan <- plan
		  .Object@replications <- replications
		  .Object@verbose <- verbose
		  .Object@parallel <- parallel
		  .Object@OK <- OK
		  .Object@only <- only
		  .Object
	  })

setMethod("getFreeVarGroup", signature("MxComputeBootstrap"),
	function(.Object) {
		result <- callNextMethod()
		for (step in c(.Object@plan)) {
			got <- getFreeVarGroup(step)
			if (length(got)) result <- append(result, got)
		}
		result
	})

setMethod("assignId", signature("MxComputeBootstrap"),
	function(.Object, id, defaultFreeSet) {
		.Object <- callNextMethod()
		defaultFreeSet <- .Object@freeSet
		id <- .Object@id
		for (sl in c('plan')) {
			slot(.Object, sl) <- assignId(slot(.Object, sl), id, defaultFreeSet)
			id <- slot(.Object, sl)@id + 1L
		}
		.Object@id <- id
		.Object
	})

setMethod("qualifyNames", signature("MxComputeBootstrap"),
	function(.Object, modelname, namespace) {
		.Object <- callNextMethod()
		for (sl in c('data')) {
			slot(.Object, sl) <- imxConvertIdentifier(slot(.Object, sl), modelname, namespace)
		}
		for (sl in c('plan')) {
			slot(.Object, sl) <- qualifyNames(slot(.Object, sl), modelname, namespace)
		}
		.Object
	})

setMethod("convertForBackend", signature("MxComputeBootstrap"),
	function(.Object, flatModel, model) {
		name <- .Object@name
		if (any(!is.integer(.Object@data))) {
			dataNum <- match(.Object@data, names(flatModel@datasets))
			if (any(is.na(dataNum))) {
				stop(paste("MxComputeBootstrap:", omxQuotes(.Object@data),
					   "not recognized as MxData"))
			}
			.Object@data <- dataNum - 1L
		}
		for (sl in c('plan')) {
			slot(.Object, sl) <- convertForBackend(slot(.Object, sl), flatModel, model)
		}
		.Object
	})

setMethod("updateFromBackend", signature("MxComputeBootstrap"),
	function(.Object, computes) {
		.Object <- callNextMethod()
		for (sl in c('plan')) {
			slot(.Object, sl) <- updateFromBackend(slot(.Object, sl), computes)
		}
		.Object
	})

mxComputeBootstrap <- function(data, plan, replications=200, ...,
			       verbose=0L, parallel=TRUE, freeSet=NA_character_,
			       OK=c("OK", "OK/green"), only=NA_integer_) {
	garbageArguments <- list(...)
	if (length(garbageArguments) > 0) {
		stop("mxComputeConfidenceInterval does not accept values for the '...' argument")
	}

	data <- vapply(data, function(e1) {
		path <- unlist(strsplit(e1, imxSeparatorChar, fixed = TRUE))
		if (length(path) == 1) {
			e1 <- paste(path, "data", sep=imxSeparatorChar)
		}
		e1
	}, "")

	new("MxComputeBootstrap", freeSet, data, plan, as.integer(replications),
	    as.integer(verbose), parallel, OK, only)
}

setMethod("displayCompute", signature(Ob="MxComputeBootstrap", indent="integer"),
	  function(Ob, indent) {
		  callNextMethod();
		  sp <- paste(rep('  ', indent), collapse="")
		  cat(sp, "$plan :", '\n')
		  displayCompute(Ob@plan, indent+1L)
		  for (sl in c("data", "replications", "OK",
			       "verbose", "parallel")) {
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

setClass(Class = "MxComputeReportExpectation",
	 contains = "BaseCompute")

setMethod("initialize", "MxComputeReportExpectation",
	  function(.Object, freeSet) {
		  .Object@name <- 'compute'
		  .Object@.persist <- TRUE
		  .Object@freeSet <- freeSet
		  .Object
	  })

##' Report expectation
##'
##' Copy the internal model expectations back to R.
##'
##' @param freeSet names of matrices containing free variables
##' @aliases
##' MxComputeReportExpectation-class

mxComputeReportExpectation <- function(freeSet=NA_character_) {
	new("MxComputeReportExpectation", freeSet)
}

#----------------------------------------------------

setClass(Class = "MxComputeSequence",
	 contains = "ComputeSteps",
	 representation = representation(
	     independent="logical"
	     ))

setMethod("initialize", "MxComputeSequence",
	  function(.Object, steps, freeSet, independent) {
		  .Object@name <- 'compute'
		  .Object@.persist <- TRUE
		  .Object@steps <- steps
		  .Object@freeSet <- freeSet
		  .Object@independent <- independent
		  .Object
	  })

##' Invoke a series of compute objects in sequence
##'
##' @param steps a list of compute objects
##' @param ... Not used; forces argument 'freeSet' to be specified by name.
##' @param freeSet Names of matrices containing free parameters.
##' @param independent Whether the steps could be executed out-of-order.
##' @aliases
##' MxComputeSequence-class
mxComputeSequence <- function(steps=list(), ..., freeSet=NA_character_, independent=FALSE) {
	garbageArguments <- list(...)
	if (length(garbageArguments) > 0) {
		stop("mxComputeSequence does not accept values for the '...' argument")
	}

	new("MxComputeSequence", steps=steps, freeSet, independent)
}

setClass(Class = "MxComputeDefault",
	 contains = "BaseCompute")

setMethod("initialize", "MxComputeDefault",
	  function(.Object, freeSet) {
		  .Object@name <- 'compute'
		  .Object@.persist <- TRUE
		  .Object@freeSet <- freeSet
		  .Object
	  })

##' Default compute plan
##'
##' The default compute plan is approximately as follows:
##' \code{mxComputeSequence(list(mxComputeGradientDescent(),
##' mxComputeConfidenceInterval(), mxComputeNumericDeriv(),
##' mxComputeStandardError(), mxComputeReportDeriv()))}
##'
##' @param freeSet names of matrices containing free variables
##' @aliases
##' MxComputeDefault-class
mxComputeDefault <- function(freeSet=NA_character_) {
	new("MxComputeDefault", freeSet)
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

setMethod("displayCompute", signature(Ob="MxComputeSequence", indent="integer"),
	  function(Ob, indent) {
		  callNextMethod();
		  sp <- paste(rep('  ', indent), collapse="")
		  cat(sp, "independent :", Ob@independent, '\n')
		  stepName <- paste("'", names(Ob@steps), "'",sep='')
		  if (length(stepName) != length(Ob@steps)) stepName <- 1:length(Ob@steps)
		  if (length(Ob@steps)) for (step in 1:length(Ob@steps)) {
			  cat(sp, "steps[[", stepName[step], "]] :", '\n')
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



omxDefaultComputePlan <- function(modelName=NULL, intervals=FALSE, useOptimizer=TRUE, 
																	optionList=options()$mxOption){
	if(length(modelName) && !is.character(modelName[1])){stop("argument 'modelName' must be a character string")}
	compute <- NULL
	fitNum <- ifelse(length(modelName), paste(modelName, 'fitfunction', sep="."), "fitfunction")
	if (!useOptimizer) {
		compute <- mxComputeSequence(list(CO=mxComputeOnce(from=fitNum, 'fit', .is.bestfit=TRUE),
																			RE=mxComputeReportExpectation()))
		} else{
		steps <- list(GD=mxComputeGradientDescent(
			fitfunction=fitNum,
			verbose=0L,	
			gradientAlgo=optionList[['Gradient algorithm']],
			gradientIterations=imxAutoOptionValue('Gradient iterations',optionList),
			gradientStepSize=imxAutoOptionValue('Gradient step size',optionList)))
			if (intervals){
				ciOpt <- mxComputeGradientDescent(
					verbose=0L,
					fitfunction=fitNum, 
					nudgeZeroStarts=FALSE,
					gradientAlgo=optionList[['Gradient algorithm']],
					gradientIterations=imxAutoOptionValue('Gradient iterations',optionList),
					gradientStepSize=imxAutoOptionValue('Gradient step size',optionList))
				cType <- ciOpt$defaultCImethod
				if (cType == 'ineq') {
					ciOpt <- mxComputeTryHard(plan=ciOpt, scale=0.05)
				}
				steps <- c(steps, CI=mxComputeConfidenceInterval(
					fitfunction=fitNum, 
					constraintType=cType,
					verbose=0L, plan=ciOpt))
			}
			if (optionList[["Calculate Hessian"]] == "Yes") {
				steps <- c(steps, ND=mxComputeNumericDeriv(
					fitfunction=fitNum, 
					stepSize=imxAutoOptionValue('Gradient step size',optionList)))
			}
			if (optionList[["Standard Errors"]] == "Yes") {
				steps <- c(steps, SE=mxComputeStandardError(), HQ=mxComputeHessianQuality())
			}
			compute <- mxComputeSequence(c(steps,
																		 RD=mxComputeReportDeriv(),
																		 RE=mxComputeReportExpectation()))
	}
	compute@.persist <- TRUE
	return(compute)
}

