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

generateThresholdColumns <- function(flatModel, model, dataName, threshName) {
	datasource <- flatModel@datasets[[dataName]]@observed
	definitionNames <- dimnames(datasource)[[2]]
	if (single.na(threshName)) {
		return(as.list(replicate(length(definitionNames), as.integer(NA))))
	}
	retval <- list()
	thresholds <- eval(substitute(mxEval(x, model, compute=TRUE),
		list(x = as.symbol(threshName))))
	threshCols <- dimnames(thresholds)[[2]]
	for(i in 1:length(definitionNames)) {
		targetName <- definitionNames[[i]]
		colIndex <- match(targetName, threshCols)
		numThresholds <- length(levels(datasource[,targetName])) - 1
		retval[[i]] <- as.integer(c(colIndex - 1, numThresholds))
	}
	return(retval)
}

verifyThresholds <- function(flatModel, model, dataName, threshName) {
	if (single.na(threshName)) {
		return()
	}
	modelName <- flatModel@name
	object <- flatModel[[threshName]]
	if (is.null(object)) {
		stop(paste("The thresholds matrix/algebra", omxQuotes(threshName), 
			"for model", omxQuotes(modelName), 
			"does not exist."), call. = FALSE)
	} else if (is(object, "MxAlgebra") || is(object, "MxMatrix")) {
	} else {
		stop(paste("The thresholds entity", omxQuotes(threshName), 
			"for model", omxQuotes(modelName), 
			"is not a MxMatrix or MxAlgebra."), call. = FALSE)
	}
	thresholds <- eval(substitute(mxEval(x, model, compute=TRUE),
		list(x = as.symbol(threshName))))
	observed <- flatModel@datasets[[dataName]]@observed
	if (is.null(dimnames(thresholds)) || is.null(dimnames(thresholds)[[2]])) {
		stop(paste("The thresholds matrix/algebra", omxQuotes(threshName), "for model", 
			omxQuotes(modelName), "does not contain column names"), call. = FALSE)
	}
	if (is.null(dimnames(observed)) || is.null(dimnames(observed)[[2]])) {
		stop(paste("The observed data matrix for model", 
			omxQuotes(modelName), "does not contain column names"), call. = FALSE)
	}
	threshNames <- dimnames(thresholds)[[2]]
	obsNames <- dimnames(observed)[[2]]
	missingNames <- setdiff(threshNames, obsNames)
	if (length(missingNames) > 0) {
		stop(paste("The following column names in the threshold",
			"matrix/algebra do not exist in the observed data matrix",
			"for model", omxQuotes(modelName), 
			":", omxQuotes(missingNames)), call. = FALSE)
	}
	for(i in 1:length(threshNames)) {
		tName <- threshNames[[i]]
		column <- thresholds[,i]
		if (!is.ordered(observed[,tName])) {
			stop(paste("In model",
				omxQuotes(modelName),
				"column",
				omxQuotes(threshNames[[i]]),
				"is not an ordered factor.",
				"Use mxFactor() on this column."), call. = FALSE)
		}
		expectedThreshCount <- length(levels(observed[,tName])) - 1
		if (nrow(thresholds) < expectedThreshCount) {
			stop(paste("In model",
				omxQuotes(modelName),
				"the number of thresholds in column",
				omxQuotes(threshNames[[i]]),
				"is less than the (l - 1), where l is equal",
				"to the number of levels in the ordinal",
				"data. Use mxFactor() on this column."), call. = FALSE)
		}
		values <- column[1:expectedThreshCount]
		sortValues <- sort(values)
		if (!all(sortValues == values)) {
			stop(paste("In model", 
				omxQuotes(modelName),
				"the thresholds in column",
				omxQuotes(tName),
				"are not in sorted order."), call. = FALSE)	
		}
	}
}

mxFactor <- function(x = character(), levels, labels = levels, exclude = NA, ordered = TRUE) {
	if(missing(levels)) {
		stop("the 'levels' argument is not optional")
	}
	if(!identical(ordered, TRUE)) {
		stop("the 'ordered' argument must be TRUE")
	}
	return(factor(x, levels, labels, exclude, ordered))
}
