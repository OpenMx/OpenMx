setClass(Class = "MxBounds",
	representation = representation(
		name = "character",
		min = "numeric",
		max = "numeric",
		parameters = "character"
	))
	
setMethod("initialize", "MxBounds",
	function(.Object, name, min, max, parameters) {
		.Object@name <- name
		.Object@min <- min
		.Object@max <- max
		.Object@parameters <- parameters
		return(.Object)
	}
)

mxBounds <- function(name = omxUntitledName(),
	min, max, parameters) {
	if (typeof(name) != "character") {
		stop(paste("Name argument is not a string",
		"(the name of the objective function)"))
	}
	if (missing(min) || !is.numeric(min)) {
		stop(paste("Min argument is not a numeric",
		"(the value of the lower bound)"))		
	}		
	if (missing(max) || !is.numeric(max)) {
		stop(paste("Max argument is not a numeric",
		"(the value of the upper bound)"))
	}
	if (missing(parameters) || 
		typeof(parameters) != "character") {
			stop(paste("Parameters argument is not a string",
		"(the vector of free parameter names)"))
	}
	if (min > max) stop(paste("min argument is greater than",
		"max argument"))
	return(new("MxBounds", name, min, max, parameters))
}