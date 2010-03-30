#
#   Copyright 2007-2010 The OpenMx Project
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

setClass(Class = "MxInterval",
	representation = representation(
		reference = "character",
		lowerdelta = "numeric",
		upperdelta = "numeric"
	))

setMethod("initialize", "MxInterval",
	function(.Object, reference, lowerdelta, upperdelta) {
		.Object@reference <- reference
		.Object@lowerdelta <- lowerdelta
		.Object@upperdelta <- upperdelta
		return(.Object)
	}
)

omxInterval <- function(reference, lowerdelta, upperdelta) {
	if (single.na(lowerdelta)) { lowerdelta <- as.numeric(NA) }
	if (single.na(upperdelta)) { upperdelta <- as.numeric(NA) }
	if (!is.character(reference) || length(reference) != 1 || is.na(reference)) {
		stop("'reference' argument must be a character string")
	}
	if (!is.numeric(lowerdelta) || length(lowerdelta) != 1) {
		stop("'lowerdelta' argument must be a numeric value")
	}
	if (!is.numeric(upperdelta) || length(upperdelta) != 1) {
		stop("'upperdelta' argument must be a numeric value")
	}
	retval <- new("MxInterval", reference, lowerdelta, upperdelta)
	return(retval)	
}

modelAddIntervals <- function(model, intervals) {
	if (length(intervals) == 0) {
		return(model)
	}
	return(model)
}

modelRemoveIntervals <- function(model, intervals) {
	if (length(intervals) == 0) {
		return(model)
	}
	return(model)
}


displayInterval <- function(object) {
	cat("MxInterval", '\n')
	cat("@reference: ", object@reference, '\n')
	cat("@lowerdelta: ", object@lowerdelta, '\n')
	cat("@upperdelta: ", object@upperdelta, '\n')
}

setMethod("print", "MxInterval", function(x,...) { displayInterval(x) })
setMethod("show", "MxInterval", function(object) { displayInterval(object) })
