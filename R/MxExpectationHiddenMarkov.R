setClass(Class = "MxExpectationHiddenMarkov",
	representation = representation(
		T0 = "MxCharOrNumber",
		G = "MxOptionalCharOrNumber",
		components = "MxOptionalCharOrNumber",
		verbose = "integer"
	),
	contains = "MxBaseExpectation")

setMethod("initialize", "MxExpectationHiddenMarkov",
          function(.Object, components, T0, G, verbose, name = 'expectation') {
		  .Object@name <- name
		  .Object@components <- components
		  .Object@T0 <- T0
		  .Object@G <- G
		  .Object@verbose <- verbose
		  .Object
	  })

setMethod("genericExpDependencies", signature("MxExpectationHiddenMarkov"),
	  function(.Object, dependencies) {
		  components <- paste(.Object@components, "expectation", sep=".")
		  sources <- c(.Object@T0, .Object@G, components)
		  dependencies <- imxAddDependency(sources, .Object@name, dependencies)
		  return(dependencies)
	  })

setMethod("genericExpFunConvert", signature("MxExpectationHiddenMarkov"), 
	  function(.Object, flatModel, model, labelsData, dependencies) {
		  if (length(.Object@components)) {
			  origComponents <- .Object@components
			  .Object@components <- vapply(.Object@components, function(group) {
				  eNumber <- match(paste(group, "expectation", sep="."),
						   names(flatModel@expectations))
				  eNumber - 1L
			  }, 1L, USE.NAMES = FALSE)
			  if (any(is.na(.Object@components))) {
				  stop(paste(name,": cannot locate expectation ",
					     omxQuotes(origComponents[is.na(.Object@components)]), sep=""),
				       call. = FALSE)
			  }
		  }
		  .Object
	  })

setMethod("qualifyNames", signature("MxExpectationHiddenMarkov"), 
	function(.Object, modelname, namespace) {
		.Object@name <- imxIdentifier(modelname, .Object@name)
		for (s in c('T0', 'G')) {
			if (is.null(slot(.Object, s))) next;
			slot(.Object, s) <-
				imxConvertIdentifier(slot(.Object, s), modelname, namespace)
		}
		.Object
})

setMethod("genericNameToNumber", signature("MxExpectationHiddenMarkov"),
	  function(.Object, flatModel, model) {
		  name <- .Object@name
		  .Object@T0 <- imxLocateIndex(flatModel, .Object@T0, name)
		  .Object@G <- imxLocateIndex(flatModel, .Object@G, name)
		  .Object
	  })

mxExpectationHiddenMarkov <- function(components, T0="T0", G=NULL, ..., verbose=0L) {
	if (length(list(...)) > 0) {
		stop(paste("Remaining parameters must be passed by name", deparse(list(...))))
	}

	new("MxExpectationHiddenMarkov", components, T0, G, as.integer(verbose))
}
