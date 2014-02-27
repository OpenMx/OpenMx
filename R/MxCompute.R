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
	     output = "list",
	     debug = "list",
	   "VIRTUAL"),
	 contains = "MxBaseNamed")

setClassUnion("MxCompute", c("NULL", "BaseCompute"))

setGeneric("convertForBackend",
	function(.Object, flatModel, model) {
		return(standardGeneric("convertForBackend"))
	})

setGeneric("updateFromBackend",
	function(.Object, computes) {
		return(standardGeneric("updateFromBackend"))
	})

setGeneric("assignId",
	function(.Object, id) {
		return(standardGeneric("assignId"))
	})

setMethod("assignId", signature("BaseCompute"),
	function(.Object, id) {
		.Object@id <- id
		.Object
	})

setGeneric("getFreeVarGroup",
	function(.Object) {
		return(standardGeneric("getFreeVarGroup"))
	})

setMethod("getFreeVarGroup", signature("BaseCompute"),
	function(.Object) {
		list()
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

#----------------------------------------------------

setClass(Class = "ComputeOperation",
	 contains = "BaseCompute",
	 representation = representation(
	   free.set = "MxOptionalChar"))

setMethod("qualifyNames", signature("ComputeOperation"),
	function(.Object, modelname, namespace) {
		.Object@name <- imxIdentifier(modelname, .Object@name)
		.Object@free.set <- imxConvertIdentifier(.Object@free.set, modelname, namespace)
		.Object
	})

setMethod("getFreeVarGroup", signature("ComputeOperation"),
	function(.Object) {
		if (length(.Object@free.set)) {
			list(.Object@id, .Object@free.set)
		} else {
			list()
		}
	})

setMethod("convertForBackend", signature("ComputeOperation"),
	function(.Object, flatModel, model) {
		name <- .Object@name
		.Object
	})

#----------------------------------------------------

setClass(Class = "MxComputeOnce",
	 contains = "ComputeOperation",
	 representation = representation(
	   what = "MxCharOrNumber",
	   verbose = "integer",
	   context = "character",
	   maxAbsChange = "logical",
	   fit = "logical",
	   gradient = "logical",
	   hessian = "logical",
	     information = "logical",
	     info.method = "MxOptionalChar",
	   ihessian = "logical",
	   hgprod = "logical",
	     .is.bestfit="logical"))

setMethod("qualifyNames", signature("MxComputeOnce"),
	function(.Object, modelname, namespace) {
		.Object <- callNextMethod()
		for (sl in c('what')) {
			slot(.Object, sl) <- imxConvertIdentifier(slot(.Object, sl), modelname, namespace)
		}
		.Object
	})

setMethod("convertForBackend", signature("MxComputeOnce"),
	function(.Object, flatModel, model) {
		.Object <- callNextMethod()
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
		   hessian, information, info.method, ihessian, hgprod, verbose, .is.bestfit) {
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
		  .Object@info.method <- info.method
		  .Object@ihessian <- ihessian
		  .Object@hgprod <- hgprod
		  .Object@.is.bestfit <- .is.bestfit
		  .Object
	  })

##' Compute something once
##'
##' The information matrix is only valid when parameters are at the
##' maximum likelihood estimate. The information matrix is returned in
##' model@output$hessian. You cannot request both the information
##' matrix and the Hessian. The information matrix is invarient to the
##' sign of the log likelihood scale whereas the Hessian is not.
##'
##' Some models are optimized for a sparse Hessian. Therefore, it can
##' be much more efficient to compute the inverse Hessian in
##' comparison to computing the Hessian and then inverting it.
##' 
##' @param what what to compute (a vector of expectation or algebra names)
##' @param free.set names of matrices containing free variables (defaults to all free variables)
##' @param context set the context for expectations
##' @param maxAbsChange compute the maximum absolute change metric
##' @param fit compute the fit statistic (often in log likelihood units)
##' @param gradient compute the analytic gradient
##' @param hessian compute the analytic Hessian
##' @param information compute the analytic information matrix
##' @param info.method which analytic method to use to compute the information
##' matrix (one of "hessian", "sandwich", "bread", and "meat")
##' @param ihessian compute the analytic inverse Hessian
##' @param hgprod not implemented
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
##'     mxComputeOnce('fitfunction', fit=TRUE))
##' factorModelFit <- mxRun(factorModel)
##' factorModelFit@output$fit  # 972.15

mxComputeOnce <- function(what, free.set=NULL, context=character(0),
			  maxAbsChange=FALSE, fit=FALSE, gradient=FALSE,
			  hessian=FALSE, information=FALSE, info.method=NULL,
			  ihessian=FALSE, hgprod=FALSE, verbose=0L, .is.bestfit=FALSE) {
	new("MxComputeOnce", what, free.set, context, maxAbsChange, fit, gradient,
	    hessian, information, info.method, ihessian, hgprod, verbose, .is.bestfit)
}

#----------------------------------------------------

setClass(Class = "MxComputeGradientDescent",
	 contains = "ComputeOperation",
	 representation = representation(
	   useGradient = "MxOptionalLogical",
	   fitfunction = "MxCharOrNumber",
	   engine = "character",
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
		.Object <- callNextMethod()
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

##' Optimize parameters using a gradient descent optimizer
##'
##' This optimizer does not require analytic derivatives of the fit
##' function. The open-source version of OpenMx only offers 1 choice,
##' CSOLNP (based on Ye, 1988).  The proprietary version of OpenMx
##' offers the choice of two optimizers, CSOLNP and NPSOL.
##'
##' @param free.set names of matrices containing free variables (defaults to all free variables)
##' @param useGradient whether to use the analytic gradient (if available)
##' @param engine specific NPSOL or CSOLNP
##' @param fitfunction name of the fitfunction (defaults to 'fitfunction')
##' @param verbose level of debugging output
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
##' factorModelFit@output$conditionNumber # 29.5

mxComputeGradientDescent <- function(free.set=NULL, useGradient=NULL,
				     engine=NULL, fitfunction='fitfunction', verbose=0L) {

	if (missing(engine)) {
		engine <- options()$mxOptions[["Default optimizer"]]
	}

	new("MxComputeGradientDescent", free.set, engine, fitfunction, useGradient, verbose)
}

#----------------------------------------------------

setClass(Class = "MxComputeNewtonRaphson",
	 contains = "ComputeOperation",
	 representation = representation(
	   fitfunction = "MxCharOrNumber",
	   maxIter = "integer",
	   tolerance = "numeric",
	   verbose = "integer",
	   carefully = "logical",
	     #output
	     iterations = "integer",
	     inform = "integer"))

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
		.Object <- callNextMethod()
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

##' Optimize parameters using the Newton-Raphson algorithm
##'
##' This optimizer requires analytic 1st and 2nd derivatives of the
##' fit function. Ramsay (1975) is used to speed convergence. Ramsay
##' can be differentially applied to different groups of parameters.
##' Comprehensive diagnostics are available by increasing the verbose
##' level. The carefully option should not be needed if the
##' derivatives are well behaved.
##'
##' @param free.set names of matrices containing free variables (defaults to all free variables)
##' @param fitfunction name of the fitfunction (defaults to 'fitfunction')
##' @param maxIter maximum number of iterations
##' @param tolerance optimization is considered converged when maximum absolute change in parameters is less than tolerance
##' @param verbose level of debugging output
##' @param carefully whether to compute the fit statistic and enforce monotonicity (defaults to FALSE)
##' @aliases
##' MxComputeNewtonRaphson-class
##' @references
##' Luenberger, D. G. & Ye, Y. (2008). \emph{Linear and nonlinear programming.} Springer.
##' 
##' Ramsay, J. O. (1975). Solving implicit equations in psychometric data analysis.
##' \emph{Psychometrika, 40}(3), 337-360.

mxComputeNewtonRaphson <- function(free.set=NULL, fitfunction='fitfunction', maxIter = 100L,
				   tolerance=1e-7, verbose=0L, carefully=FALSE)
{
	new("MxComputeNewtonRaphson", free.set, fitfunction, maxIter, tolerance, verbose, carefully)
}

#----------------------------------------------------

setClass(Class = "ComputeSteps",
	 contains = "BaseCompute",
	 representation = representation(
	   steps = "list"))

setMethod("getFreeVarGroup", signature("ComputeSteps"),
	function(.Object) {
		result <- list()
		for (step in .Object@steps) {
			got <- getFreeVarGroup(step)
			if (length(got)) result <- append(result, got)
		}
		result
	})

setMethod("assignId", signature("ComputeSteps"),
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

setMethod("qualifyNames", signature("ComputeSteps"),
	function(.Object, modelname, namespace) {
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
	  function(.Object, steps, maxIter, tolerance, verbose) {
		  .Object@name <- 'compute'
		  .Object@steps <- steps
		  .Object@maxIter <- maxIter
		  .Object@tolerance <- tolerance
		  .Object@verbose <- verbose
		  .Object
	  })

##' Repeatedly invoke a series of compute objects until change is less than tolerance
##'
##' One step (typically the last) must compute the fit or maxAbsChange.
##'
##' @param steps a list of compute objects
##' @param maxIter the maximum number of iterations
##' @param tolerance iterates until change is less than tolerance
##' @param verbose level of debugging output
##' @aliases
##' MxComputeIterate-class
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
	 contains = "ComputeOperation",
	 representation = representation(
	     what = "MxCharOrNumber",
	     mstep.fit = "MxCompute",
	     fit = "MxCompute",
	     maxIter = "integer",
	     tolerance = "numeric",
	     verbose = "integer",
	     ramsay="logical",
	     information="logical",
	     semMethod="MxOptionalNumeric",
	     semDebug="logical",
	     noiseTarget="numeric",
	     noiseTolerance="numeric",
	     info.method="character",
	     semFixSymmetry="logical",
	     agileMaxIter="integer"))

setMethod("assignId", signature("MxComputeEM"),
	function(.Object, id) {
		.Object@mstep.fit <- assignId(.Object@mstep.fit, id)
		.Object@fit <- assignId(.Object@fit, id + 1L)
		.Object@id <- id + 2L
		.Object
	})

setMethod("getFreeVarGroup", signature("MxComputeEM"),
	function(.Object) {
		result <- callNextMethod()
		for (step in c(.Object@mstep.fit, .Object@fit)) {
			got <- getFreeVarGroup(step)
			if (length(got)) result <- append(result, got)
		}
		result
	})

setMethod("qualifyNames", signature("MxComputeEM"),
	function(.Object, modelname, namespace) {
		.Object <- callNextMethod()
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
		.Object <- callNextMethod()
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

setMethod("updateFromBackend", signature("MxComputeEM"),
	function(.Object, computes) {
		.Object <- callNextMethod()
		.Object@mstep.fit <- updateFromBackend(.Object@mstep.fit, computes)
		.Object@fit <- updateFromBackend(.Object@fit, computes)
		.Object
	})

setMethod("initialize", "MxComputeEM",
	  function(.Object, what, mstep.fit, fit, maxIter, tolerance,
		   verbose, ramsay, information, noiseTarget, noiseTolerance, semDebug,
		   semMethod, info.method, semFixSymmetry, agileMaxIter) {
		  .Object@name <- 'compute'
		  .Object@what <- what
		  .Object@mstep.fit <- mstep.fit
		  .Object@fit <- fit
		  .Object@maxIter <- maxIter
		  .Object@tolerance <- tolerance
		  .Object@verbose <- verbose
		  .Object@ramsay <- ramsay
		  .Object@information <- information
		  .Object@noiseTarget <- noiseTarget
		  .Object@noiseTolerance <- noiseTolerance
		  .Object@semDebug <- semDebug
		  .Object@semMethod <- semMethod
		  .Object@info.method <- info.method
		  .Object@semFixSymmetry <- semFixSymmetry
		  .Object@agileMaxIter <- agileMaxIter
		  .Object
	  })

mxComputeEM <- function(expectation, mstep.fit, fit, maxIter=500L, tolerance=1e-4,
			verbose=0L, ramsay=TRUE, information=FALSE, noiseTarget=exp(-3), noiseTolerance=exp(1.5),
			semDebug=FALSE, semMethod=NULL, info.method="hessian", semFixSymmetry=TRUE, agileMaxIter=1L) {
	new("MxComputeEM", what=expectation, mstep.fit, fit, maxIter=maxIter,
	    tolerance=tolerance, verbose, ramsay, information, noiseTarget,
	    noiseTolerance, semDebug, semMethod, info.method, semFixSymmetry, agileMaxIter)
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

setClass(Class = "MxComputeNumericDeriv",
	 contains = "ComputeOperation",
	 representation = representation(
	   fitfunction = "MxCharOrNumber",
	     parallel = "logical",
	     stepSize = "numeric",
	     iterations = "integer"))

setMethod("qualifyNames", signature("MxComputeNumericDeriv"),
	function(.Object, modelname, namespace) {
		.Object <- callNextMethod();
		.Object@fitfunction <- imxConvertIdentifier(.Object@fitfunction, modelname, namespace)
		.Object
	})

setMethod("convertForBackend", signature("MxComputeNumericDeriv"),
	function(.Object, flatModel, model) {
		.Object <- callNextMethod();
		name <- .Object@name
		if (is.character(.Object@fitfunction)) {
			.Object@fitfunction <- imxLocateIndex(flatModel, .Object@fitfunction, name)
		}
		.Object
	})

setMethod("initialize", "MxComputeNumericDeriv",
	  function(.Object, free.set, fit, parallel, stepSize, iterations) {
		  .Object@name <- 'compute'
		  .Object@free.set <- free.set
		  .Object@fitfunction <- fit
		  .Object@parallel <- parallel
		  .Object@stepSize <- stepSize
		  .Object@iterations <- iterations
		  .Object
	  })

##' Numerically estimate Hessian using Richardson extrapolation
##'
##' For N free parameters, Richardson extrapolation requires
##' (iterations * (N^2 + N)) function evaluations.
##' 
##' The implementation is closely based on the numDeriv R package.
##' 
##' @param free.set names of matrices containing free variables (defaults to all free variables)
##' @param fitfunction name of the fitfunction (defaults to 'fitfunction')
##' @param parallel whether to evaluate the fitfunction in parallel (defaults to TRUE)
##' @param stepSize starting set size (defaults to 0.0001)
##' @param iterations number of Richardson extrapolation iterations (defaults to 4L)
##' @aliases
##' MxComputeNumericDeriv-class
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
##'     mxComputeNumericDeriv())
##' factorModelFit <- mxRun(factorModel)
##' factorModelFit@output$hessian

mxComputeNumericDeriv <- function(free.set=NULL, fitfunction='fitfunction',
				      parallel=TRUE, stepSize=0.0001, iterations=4L) {
	new("MxComputeNumericDeriv", free.set, fitfunction, parallel, stepSize, iterations)
}

#----------------------------------------------------

setClass(Class = "MxComputeStandardError",
	 representation = representation(
	     forcePositiveDefinite = "logical"),
	 contains = "ComputeOperation")

setMethod("initialize", "MxComputeStandardError",
	  function(.Object, free.set, forcePositiveDefinite) {
		  .Object@name <- 'compute'
		  .Object@free.set <- free.set
		  .Object@forcePositiveDefinite <- forcePositiveDefinite
		  .Object
	  })

##' Compute standard errors given the Hessian or inverse Hessian
##'
##' @param free.set names of matrices containing free variables (defaults to all free variables)
##' @aliases
##' MxComputeStandardError-class

mxComputeStandardError <- function(free.set=NULL, forcePositiveDefinite=TRUE) {
	new("MxComputeStandardError", free.set, forcePositiveDefinite)
}

#----------------------------------------------------

setClass(Class = "MxComputeHessianQuality",
	 contains = "ComputeOperation")

setMethod("initialize", "MxComputeHessianQuality",
	  function(.Object, free.set) {
		  .Object@name <- 'compute'
		  .Object@free.set <- free.set
		  .Object
	  })

##' Compute the quality of the Hessian
##'
##' Tests whether the Hessian is positive definite
##' (model@output$infoDefinite) and, if so, computes the condition
##' number (model@output$conditionNumber). See Luenberger & Ye (2008)
##' Second Order Test (p. 190) and Condition Number (p. 239).
##' 
##' @param free.set names of matrices containing free variables (defaults to all free variables)
##' @aliases
##' MxComputeHessianQuality-class
##' @references
##' Luenberger, D. G. & Ye, Y. (2008). Linear and nonlinear programming. Springer.

mxComputeHessianQuality <- function(free.set=NULL) {
	new("MxComputeHessianQuality", free.set)
}

#----------------------------------------------------

setClass(Class = "MxComputeSequence",
	 contains = "ComputeSteps")

setMethod("initialize", "MxComputeSequence",
	  function(.Object, steps) {
		  .Object@name <- 'compute'
		  .Object@steps <- steps
		  .Object
	  })

##' Invoke a series of compute objects in sequence
##'
##' @param steps a list of compute objects
##' @aliases
##' MxComputeSequence-class
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

displayComputeOperation <- function(opt) {
	cat(class(opt), omxQuotes(opt@name), '\n')
	cat("@id :", opt@id, '\n')
	cat("@free.set :", omxQuotes(opt@free.set), '\n')
	invisible(opt)
}

setMethod("print", "ComputeOperation", function(x, ...) displayComputeOperation(x))
setMethod("show",  "ComputeOperation", function(object) displayComputeOperation(object))

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
	if (is.null(flatModel@compute)) return()
	convertForBackend(flatModel@compute, flatModel, model)
}

updateModelCompute <- function(model, computes) {
	if (is.null(model@compute)) return()
	updateFromBackend(model@compute, computes)
}
