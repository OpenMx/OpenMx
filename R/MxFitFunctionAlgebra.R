#
#   Copyright 2007-2015 The OpenMx Project
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


setClass(Class = "MxFitFunctionAlgebra",
	representation = representation(
		algebra = "MxCharOrNumber",
	    units = "character",
		numObs = "numeric",
		numStats = "numeric",
	    gradient = "MxCharOrNumber",
	    hessian = "MxCharOrNumber",
	    verbose = "integer"),
	contains = "MxBaseFitFunction")

setMethod("initialize", "MxFitFunctionAlgebra",
	function(.Object, algebra, units, numObs, numStats, gradient, hessian, verbose, name = 'fitfunction') {
		.Object@name <- name
		.Object@algebra <- algebra
		.Object@units <- units
		.Object@numObs <- numObs
		.Object@numStats <- numStats
		.Object@gradient <- gradient
		.Object@hessian <- hessian
		.Object@verbose <- verbose
		return(.Object)
	}
)

setMethod("genericFitDependencies", signature("MxFitFunctionAlgebra"),
	function(.Object, flatModel, dependencies) {
	dependencies <- callNextMethod()
	for (sl in c('algebra', 'gradient', 'hessian')) {
		thing <- slot(.Object, sl)
		if (is.na(thing)) next
		dependencies <- imxAddDependency(thing, .Object@name, dependencies)
	}
	return(dependencies)
})

setMethod("genericFitFunConvert", signature("MxFitFunctionAlgebra"), 
	function(.Object, flatModel, model, labelsData, defVars, dependencies) {
		name <- .Object@name
		algebra <- .Object@algebra
		if (is.na(algebra) && is.na(.Object@gradient) && is.na(.Object@hessian)) {
			modelname <- imxReverseIdentifier(model, .Object@name)[[1]]
			msg <- paste("The algebra name cannot be NA",
			"for the algebra fit function of model", omxQuotes(modelname))
			stop(msg, call. = FALSE)
		}
		modelname <- imxReverseIdentifier(model, .Object@name)[[1]]
		expectName <- paste(modelname, "expectation", sep=".")
		if (expectName %in% names(flatModel@expectations)) {
			expectIndex <- imxLocateIndex(flatModel, expectName, name)
		} else {
			expectIndex <- as.integer(NA)
		}
		.Object@expectation <- expectIndex
		for (sl in c('algebra', 'gradient', 'hessian')) {
			slot(.Object, sl) <- imxLocateIndex(flatModel, slot(.Object, sl), name)
		}
		return(.Object)
})

setMethod("qualifyNames", signature("MxFitFunctionAlgebra"), 
	function(.Object, modelname, namespace) {
		.Object@name <- imxIdentifier(modelname, .Object@name)
		for (sl in c('algebra', 'gradient', 'hessian')) {
			slot(.Object, sl) <- imxConvertIdentifier(slot(.Object, sl), modelname, namespace)
		}
		return(.Object)
})

setMethod("genericFitRename", signature("MxFitFunctionAlgebra"),
	function(.Object, oldname, newname) {
		for (sl in c('algebra', 'gradient', 'hessian')) {
			slot(.Object, sl) <- renameReference(slot(.Object, sl), oldname, newname)
		}
		return(.Object)
})

setMethod("generateReferenceModels", "MxFitFunctionAlgebra",
	function(.Object, model) {
		msg <- paste("Don't know how to make reference models for a model with a ",
			     class(.Object), " fit function.", sep="")
		msg <- paste(msg, "\n",
			     "If you're using this for a mutligroup model, very likely, you can replace your mxFitFunctionAlgebra() call with", "\n",
			     "mxFitFunctionMultigroup(c('submodelName1', 'submodelName2', ...))", "\n\n",
			     "See ?mxFitFunctionMultigroup() to learn more.", sep="")
		stop(msg)
	})

mxFitFunctionAlgebra <- function(algebra, numObs = NA, numStats = NA, ...,
				 gradient=NA_character_, hessian=NA_character_,
				 verbose=0L, units="-2lnL")
{
	garbageArguments <- list(...)
	if (length(garbageArguments) > 0) {
		stop("mxFitFunctionAlgebra does not accept values for the '...' argument")
	}

	if (is.null(algebra)) {
		algebra <- NA_character_
	} else if (missing(algebra) || typeof(algebra) != "character") {
		stop("Algebra argument is not a string (the name of the algebra)")
	}
	if (single.na(numObs)) {
		numObs <- as.numeric(NA)
	}
	if (single.na(numStats)) {
		numStats <- as.numeric(NA)
	}
	return(new("MxFitFunctionAlgebra", algebra, units, numObs, numStats, gradient, hessian, verbose))
}

displayMxFitFunctionAlgebra <- function(fitfunction) {
	cat("MxFitFunctionAlgebra", omxQuotes(fitfunction@name), '\n')
	cat("$algebra: ", omxQuotes(fitfunction@algebra), '\n')	
	cat("$units: ", omxQuotes(fitfunction@units), '\n')	
	cat("$numObs: ", fitfunction@numObs, '\n')
	cat("$numStats: ", fitfunction@numStats, '\n')
	if (length(fitfunction@result) == 0) {
		cat("$result: (not yet computed) ")
	} else {
		cat("$result:\n")
	}
	print(fitfunction@result)
	invisible(fitfunction)
}

setMethod("print", "MxFitFunctionAlgebra", function(x,...) { 
	displayMxFitFunctionAlgebra(x) 
})

setMethod("show", "MxFitFunctionAlgebra", function(object) { 
	displayMxFitFunctionAlgebra(object) 
})
