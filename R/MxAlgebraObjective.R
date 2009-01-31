setClass(Class = "MxAlgebraObjective",
	representation = representation(
		algebra = "MxCharOrNumber"),
	contains = "MxObjective")

setMethod("initialize", "MxAlgebraObjective",
	function(.Object, name, algebra) {
		.Object@name <- name
		.Object@algebra <- algebra
		return(.Object)
	}
)

setMethod("omxObjFunConvert", signature("MxAlgebraObjective", "MxModel"), function(.Object, model) {
		name <- .Object@name
		algebra <- .Object@algebra
		algebraIndex <- omxLocateIndex(model, algebra)
		if (is.na(algebraIndex)) {
			stop(paste("Could not find a matrix/algebra with name", algebra, "in the model."))
		}
		return(new("MxAlgebraObjective", name, algebraIndex))
})


mxAlgebraObjective <- function(name = omxUntitledName(), algebra) {
	if (typeof(name) != "character") {
		stop("Name argument is not a string (the name of the objective function)")
	}
	if (typeof(algebra) != "character") {
		stop("Algebra argument is not a string (the name of the algebra)")
	}
	return(new("MxAlgebraObjective", name, algebra))
}