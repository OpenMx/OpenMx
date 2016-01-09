#
#   Copyright 2007-2016 The OpenMx Project
#
#   Licensed under the Apache License, Version 2.0 (the "License");
#   you may not use this file except in compliance with the License.
#   You may obtain a copy of the License at
# 
#        http://www.apache.org/licenses/LICENSE-2.0
# 
#   Unless required by applicable law or agreed to in writing, software
#   distributed under the License is distributed on an "AS IS" BASIS,
#   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#   See the License for the specific language governing permissions and
#   limitations under the License.


setClass(Class = "MxBounds",
	representation = representation(
		min = "numeric",
		max = "numeric",
		parameters = "character"
	))
	
setMethod("initialize", "MxBounds",
	function(.Object, min, max, parameters) {
		.Object@min <- min
		.Object@max <- max
		.Object@parameters <- parameters
		return(.Object)
	}
)

mxBounds <- function(parameters, min = NA, max = NA) {
	if (length(min) > 1) {
		stop("Only a single min value may be specified")
	}
	if (length(max) > 1) {
		stop("Only a single max value may be specified")
	}
	if(is.na(min)) min <- NA_real_
	if(is.na(max)) max <- NA_real_
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
	return(new("MxBounds", min, max, parameters))
}

checkBounds <- function(model, bounds) {
	parameters <- omxGetParameters(model)
	parameters <- names(parameters)
	for(i in 1:length(bounds)) {
		bound <- bounds[[i]]
		nomatch <- !(bound@parameters %in% parameters)
		if (any(nomatch)) {
			stop(paste("In model ", omxQuotes(model@name),
				" the following parameter names are",
				" used in mxBounds() but do",
				" not exist in the model: ",
				omxQuotes(bound@parameters[nomatch]), '.',
				sep = ''), call. = FALSE)
		}
	}
}

modelAddBounds <- function(model, bounds) {
	if (length(bounds) == 0) {
		return(model)
	}
	checkBounds(model, bounds)
	return(modelAddBoundsHelper(model, bounds))
}

modelAddBoundsHelper <- function(model, bounds) {
	if (length(model@matrices) > 0) {
		for(i in 1:length(bounds)) {
			for(j in 1:length(model@matrices)) {
				matrix <- model@matrices[[j]]
				matches <- matrix$labels %in% bounds[[i]]@parameters
				if (any(matches)) {					
					matrix$lbound[matches] <- bounds[[i]]@min
					matrix$ubound[matches] <- bounds[[i]]@max
					model@matrices[[j]] <- matrix
				}
			}
		}
	}
	if (length(model@submodels) > 0) {
		for(i in 1:length(model@submodels)) {
			model@submodels[[i]] <- modelAddBoundsHelper(model@submodels[[i]], bounds)
		}
	}
	return(model)
}

modelRemoveBounds <- function(model, bounds) {
	if (length(bounds) == 0) {
		return(model)
	}
	if (length(model@matrices) > 0) {
		for(i in 1:length(bounds)) {
			for(j in 1:length(model@matrices)) {
				matrix <- model@matrices[[j]]
				matches <- matrix@labels %in% bounds[[i]]@parameters
				matrix@lbound[matches] <- as.numeric(NA)
				matrix@ubound[matches] <- as.numeric(NA)
				model@matrices[[j]] <- matrix
			}
		}
	}
	if (length(model@submodels) > 0) {
		for(i in 1:length(model@submodels)) {
			model@submodels[[i]] <- modelRemoveBounds(model@submodels[[i]], bounds)
		}
	}
	return(model)
}
