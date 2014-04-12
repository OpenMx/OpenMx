setClass(Class = "MxFitFunctionMultigroup",
	 representation = representation(
	   groups = "MxOptionalCharOrNumber",
	     verbose= "integer"),
	 contains = "MxBaseFitFunction")

setMethod("initialize", "MxFitFunctionMultigroup",
	function(.Object, groups, verbose, name = 'fitfunction') {
		.Object@name <- name
		.Object@groups <- groups
		.Object@verbose <- verbose
		return(.Object)
	}
)

setMethod("genericFitDependencies", signature("MxFitFunctionMultigroup"),
	function(.Object, flatModel, dependencies) {
	dependencies <- callNextMethod()
	dependencies <- imxAddDependency(.Object@groups, .Object@name, dependencies)
	return(dependencies)
})

setMethod("qualifyNames", signature("MxFitFunctionMultigroup"), 
	function(.Object, modelname, namespace) {
		.Object@name <- imxIdentifier(modelname, .Object@name)
		return(.Object)
})

# "model.algebra" or "model" for "model.fitfunction"
setMethod("genericFitFunConvert", "MxFitFunctionMultigroup", 
	function(.Object, flatModel, model, labelsData, defVars, dependencies) {
		name <- .Object@name
		if (length(.Object@groups)) .Object@groups <- vapply(.Object@groups, function(group) {
			path <- unlist(strsplit(group, imxSeparatorChar, fixed = TRUE))
			if (length(path) == 1) {
				ff <- paste(path, "fitfunction", sep=".")
				length(model@algebras) + imxLocateIndex(flatModel, ff, name)
			} else if (length(path) == 2) {
				# restrict to algebra or fitfunction TODO
				imxLocateIndex(flatModel, group, name)
			}
		}, 1L)
		return(.Object)
})

##' Aggregate fit statistics from submodels
##'
##' This is a simple fit function that sums the fit statistics
##' from other fit functions, typically in submodels. It is almost
##' equivalent to,
##'
##' \code{mxAlgebra(model1.objective + model2.objective, name="alg")}
##'
##' and
##'
##' \code{mxFitFunctionAlgebra("alg")}
##'
##' However, in addition to the fit statistic, mxFitFunctionMultigroup
##' also aggregates analytic derivative calculations.
##'
##' @param groups vector of fit function names
##' @param ...  Not used.  Forces remaining arguments to be specified by name.
##' @param verbose the level of debugging output
##' @aliases
##' MxFitFunctionMultigroup-class
mxFitFunctionMultigroup <- function(groups, ..., verbose=0L) {
	garbageArguments <- list(...)
	if (length(garbageArguments) > 0) {
		stop("mxFitFunctionMultigroup does not accept values for the '...' argument")
	}

	if (length(groups) == 0) stop("mxFitFunctionMultigroup: at least 1 fitfunction must be provided")

	return(new("MxFitFunctionMultigroup", groups, verbose))
}
