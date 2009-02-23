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
	min = NA, max = NA, parameters) {
	if(is.na(min)) min <- NA_real_
	if(is.na(max)) max <- NA_real_
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
	if (!is.na(min) && !is.na(max) && min > max) { 
		msg <- paste("min argument is greater than",
			"max argument")
		stop(msg)
	}
	return(new("MxBounds", name, min, max, parameters))
}

omxLocateBounds <- function(bounds, parameterName) {
	filter <- lapply(bounds, function(x) {
		if (parameterName %in% x@parameters) {
			return(c(x@min, x@max))
		} else {return(NA)}
	})
	filter <- filter[!is.na(filter)]
	if (length(filter) == 0) {
		return(c(NA_real_,NA_real_))
	} else if (length(filter) == 1) {
		return(filter[[1]])
	} else {
		msg <- paste("The parameter", omxQuotes(parameterName),
			"has multiple specifications of bounds in the model")
		stop(msg, call.=FALSE)
	}
}