setClass(Class = "MxExpectationMixture",
	representation = representation(
		weights = "MxCharOrNumber",
		components = "MxOptionalCharOrNumber",
		verbose = "integer",
		scale = "character"
	),
	contains = "MxBaseExpectation")

setMethod("initialize", "MxExpectationMixture",
          function(.Object, components, weights, verbose, scale, name = 'expectation') {
		  .Object@data = as.integer(NA)
		  .Object@name <- name
		  .Object@components <- components
		  .Object@weights <- weights
		  .Object@verbose <- verbose
		  .Object@scale <- scale
		  .Object
	  })

setMethod("genericExpDependencies", signature("MxExpectationMixture"),
	  function(.Object, dependencies) {
		  components <- paste(.Object@components, "expectation", sep=".")
		  sources <- c(.Object@weights, components)
		  dependencies <- imxAddDependency(sources, .Object@name, dependencies)
		  return(dependencies)
	  })

setMethod("genericExpFunConvert", signature("MxExpectationMixture"), 
	  function(.Object, flatModel, model, labelsData, dependencies) {
		  if (length(.Object@components)) {
			  origComponents <- .Object@components
			  .Object@components <- vapply(.Object@components, function(group) {
				  eNumber <- match(paste(group, "expectation", sep="."),
						   names(flatModel@expectations))
				  eNumber - 1L
			  }, 1L, USE.NAMES = FALSE)
			  if (any(is.na(.Object@components))) {
				  stop(paste(model@name,": cannot locate expectation ",
					     omxQuotes(origComponents[is.na(.Object@components)]), sep=""),
				       call. = FALSE)
			  }
		  }
		  .Object
	  })

setMethod("qualifyNames", signature("MxExpectationMixture"), 
	function(.Object, modelname, namespace) {
		.Object@name <- imxIdentifier(modelname, .Object@name)
		.Object@data <- imxConvertIdentifier(.Object@data, 
			modelname, namespace)
		for (s in c('weights')) {
			if (is.null(slot(.Object, s))) next;
			slot(.Object, s) <-
				imxConvertIdentifier(slot(.Object, s), modelname, namespace)
		}
		.Object
})

setMethod("genericNameToNumber", signature("MxExpectationMixture"),
	  function(.Object, flatModel, model) {
		  name <- .Object@name
		  .Object@data <- imxLocateIndex(flatModel, .Object@data, name)
		  .Object@weights <- imxLocateIndex(flatModel, .Object@weights, name)
		  .Object
	  })

mxExpectationMixture <- function(components, weights="weights",
				      ..., verbose=0L, scale=c('softmax','sum')) {
	if (length(list(...)) > 0) {
		stop(paste("Remaining parameters must be passed by name", deparse(list(...))))
	}

	scale <- match.arg(scale)

	new("MxExpectationMixture", components, weights,
	    as.integer(verbose), scale)
}
