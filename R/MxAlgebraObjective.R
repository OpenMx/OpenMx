setClass(Class = "MxAlgebraObjective",
	representation = representation(
		algebra = "numeric"),
	contains = "MxObjective")

setMethod("initialize", "MxAlgebraObjective",
	function(.Object, algebra) {
		.Object@algebra <- algebra
		return(.Object)
	}
)


mxAlgebraObjective <- function(model, algebra) {
	if(class(model)[[1]] != "MxModel") {
		stop("First argument is not an MxModel object")
	}
	if (typeof(algebra) != "character") {
		stop("Second argument is not a string (the name of the algebra)")
	}
	algebraIndex <- omxLocateIndex(model, algebra)
	if (is.na(algebraIndex)) {
		stop(paste("Could not find a matrix/algebra with name", algebra, "in the model."))
	}
	return(new("MxAlgebraObjective", algebraIndex))
}