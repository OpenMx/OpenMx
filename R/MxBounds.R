#
#   Copyright 2007-2009 The OpenMx Project
#
#   Licensed under the Apache License, Version 2.0 (the "License");
#   you may not use this file except in compliance with the License.
#   You may obtain a copy of the License at
# 
#        http://www.apache.org/licenses/LICENSE-2.0
# 
#    Unless required by applicable law or agreed to in writing, software
#    distributed under the License is distributed on an "AS IS" BASIS,
#    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#   See the License for the specific language governing permissions and
#   limitations under the License.


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

mxBounds <- function(parameters, min = NA, max = NA, name = NA) {
	if(is.na(min)) min <- NA_real_
	if(is.na(max)) max <- NA_real_
	if(is.na(name)) name <- omxUntitledName()
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
