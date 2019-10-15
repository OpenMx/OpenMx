setClass(Class = "MxExpectationHiddenMarkov",
	representation = representation(
		initial = "MxCharOrNumber",
		transition = "MxOptionalCharOrNumber",
		components = "MxOptionalCharOrNumber",
		verbose = "integer",
		scale = "character"
	),
	contains = "MxBaseExpectation")

setMethod("initialize", "MxExpectationHiddenMarkov",
          function(.Object, components, initial, transition, verbose, scale, name = 'expectation') {
		  .Object@data = as.integer(NA)
		  .Object@name <- name
		  .Object@components <- components
		  .Object@initial <- initial
		  .Object@transition <- transition
		  .Object@verbose <- verbose
		  .Object@scale <- scale
		  .Object
	  })

setMethod("genericExpDependencies", signature("MxExpectationHiddenMarkov"),
	  function(.Object, dependencies) {
		  components <- paste(.Object@components, "expectation", sep=".")
		  sources <- c(.Object@initial, .Object@transition, components)
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
				  stop(paste(model@name,": cannot locate expectation ",
					     omxQuotes(origComponents[is.na(.Object@components)]), sep=""),
				       call. = FALSE)
			  }
		  }
		  .Object
	  })

setMethod("qualifyNames", signature("MxExpectationHiddenMarkov"), 
	function(.Object, modelname, namespace) {
		.Object@name <- imxIdentifier(modelname, .Object@name)
		.Object@data <- imxConvertIdentifier(.Object@data, 
			modelname, namespace)
		for (s in c('initial', 'transition')) {
			if (is.null(slot(.Object, s))) next;
			slot(.Object, s) <-
				imxConvertIdentifier(slot(.Object, s), modelname, namespace)
		}
		.Object
})

setMethod("genericNameToNumber", signature("MxExpectationHiddenMarkov"),
	  function(.Object, flatModel, model) {
		  name <- .Object@name
		  .Object@data <- imxLocateIndex(flatModel, .Object@data, name)
		  .Object@initial <- imxLocateIndex(flatModel, .Object@initial, name)
		  .Object@transition <- imxLocateIndex(flatModel, .Object@transition, name)
		  .Object
	  })

mxExpectationHiddenMarkov <- function(components, initial="initial", transition=NULL,
				      ..., verbose=0L, scale=c('softmax', 'sum', 'none')) {
  prohibitDotdotdot(list(...))

	scale <- match.arg(scale)

	new("MxExpectationHiddenMarkov", components, initial, transition,
	    as.integer(verbose), scale)
}
