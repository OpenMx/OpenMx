setClass(Class = "MxAlgebraObjective",
	representation = representation(
		algebra = "MxCharOrNumber"),
	contains = "MxBaseObjective")

setMethod("initialize", "MxAlgebraObjective",
	function(.Object, name, algebra) {
		.Object@name <- name
		.Object@algebra <- algebra
		return(.Object)
	}
)

setMethod("omxObjFunConvert", signature("MxAlgebraObjective", "MxFlatModel"), function(.Object, model) {
		name <- .Object@name
		algebra <- .Object@algebra
		algebraIndex <- omxLocateIndex(model, algebra, name)
		if (is.na(algebraIndex)) {
			stop(paste("Could not find a matrix/algebra with name", 
				algebra, "in the model."))
		}
		return(new("MxAlgebraObjective", name, algebraIndex))
})


mxAlgebraObjective <- function(algebra, name = omxUntitledName()) {
	if (typeof(name) != "character") {
		stop("Name argument is not a string (the name of the objective function)")
	}
	if (missing(algebra) || typeof(algebra) != "character") {
		stop("Algebra argument is not a string (the name of the algebra)")
	}
	return(new("MxAlgebraObjective", name, algebra))
}
