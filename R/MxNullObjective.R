setClass(Class = "MxNullObjective",
	representation = representation(),
	contains = "MxObjective")

setMethod("initialize", "MxNullObjective",
	function(.Object, name) {
		.Object@name <- name
		return(.Object)
	}
)

setMethod("omxObjFunConvert", signature("MxNullObjective", "MxModel"), function(.Object, model) {
		return(.Object)
	}
)

mxNullObjective <- function(name = omxUntitledName()) {
	return(new("MxNullObjective", name))
}