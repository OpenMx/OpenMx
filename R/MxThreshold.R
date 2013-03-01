#
#   Copyright 2007-2013 The OpenMx Project
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

generateThresholdColumns <- function(flatModel, model, labelsData, covarianceColumnNames, dataName, threshName) {
	covarianceLength <- length(covarianceColumnNames)
	thresholdColumns <- replicate(covarianceLength, as.integer(NA))
	thresholdLevels <- replicate(covarianceLength, as.integer(NA))
	if (single.na(threshName)) {
		return(list(thresholdColumns, thresholdLevels))
	}

	tuple <- evaluateMxObject(threshName, flatModel, labelsData, new.env(parent = emptyenv()))
	thresholds <- tuple[[1]]

	thresholdColumnNames <- dimnames(thresholds)[[2]]
	datasource <- flatModel@datasets[[dataName]]@observed
	for(i in 1:length(thresholdColumnNames)) {
		oneThresholdName <- thresholdColumnNames[[i]]
		numThresholds <- length(levels(datasource[, oneThresholdName])) - 1
		columnIndex <- match(oneThresholdName, covarianceColumnNames)
		thresholdColumns[[columnIndex]] <- as.integer(i - 1)
		thresholdLevels[[columnIndex]] <- as.integer(numThresholds)
	}
	return(list(thresholdColumns, thresholdLevels))
}

verifyThresholds <- function(flatModel, model, labelsData, dataName, covNames, threshName) {
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
	tuple <- evaluateMxObject(threshName, flatModel, labelsData, new.env(parent = emptyenv()))
	thresholds <- tuple[[1]]
	observed <- flatModel@datasets[[dataName]]@observed
	if (is.null(dimnames(thresholds)) || is.null(dimnames(thresholds)[[2]])) {
		stop(paste("The thresholds matrix/algebra", omxQuotes(threshName), "for model", 
			omxQuotes(modelName), "does not contain column names"), call. = FALSE)
	}
	if (!is.data.frame(observed)) {
		stop(paste("The observed data for model", 
			omxQuotes(modelName), "is not a data.frame object"), call. = FALSE)
	}
	if (is.null(dimnames(observed)) || is.null(dimnames(observed)[[2]])) {
		stop(paste("The observed data frame for model", 
			omxQuotes(modelName), "does not contain column names"), call. = FALSE)
	}
	threshNames <- dimnames(thresholds)[[2]]
	obsNames <- dimnames(observed)[[2]]
	missingNames <- setdiff(threshNames, obsNames)
	if (length(missingNames) > 0) {
		stop(paste("The column name(s)", omxQuotes(missingNames), 
			"appear in the thresholds matrix",
			omxQuotes(simplifyName(threshName, modelName)), "but not",
			"in the observed data.frame or matrix",
			"in model", omxQuotes(modelName)), call. = FALSE)
	}
	missingNames <- setdiff(threshNames, covNames)
	if (length(missingNames) > 0) {
		stop(paste("The column name(s)", omxQuotes(missingNames), 
			"appear in the thresholds matrix",
			omxQuotes(simplifyName(threshName, modelName)), "but not",
			"the manifest variables of model", omxQuotes(modelName)), call. = FALSE)
	}
	for(i in 1:length(threshNames)) {
		tName <- threshNames[[i]]
		tColumn <- thresholds[,i]
		observedColumn <- observed[,tName]
		if (!all(is.na(observedColumn))) {
			if (!is.ordered(observedColumn)) {
				stop(paste("In model",
					omxQuotes(modelName),
					"column",
					omxQuotes(tName),
					"is not an ordered factor.",
					"Use mxFactor() on this column."), call. = FALSE)
			}
			expectedThreshCount <- length(levels(observedColumn)) - 1
			if (nrow(thresholds) < expectedThreshCount) {
				stop(paste("In model",
					omxQuotes(modelName),
					"the number of thresholds in column",
					omxQuotes(tName),
					"is less than the (l - 1), where l is equal",
					"to the number of levels in the ordinal",
					"data. Use mxFactor() on this column."), call. = FALSE)
			}
			values <- tColumn[1:expectedThreshCount]
			if (any(is.na(values))) {
				stop(paste("In model", 
					omxQuotes(modelName),
					"I was expecting", expectedThreshCount, 
					"thresholds in column", omxQuotes(tName), 
					"of matrix/algebra", 
					omxQuotes(simplifyName(threshName, modelName)),
					"but I hit NA values after only", 
					length(values[!is.na(values)]), "thresholds.",
					"You need to increase the number of",
					"free thresholds for", omxQuotes(tName), 
					"and give them values other than NA"), call. = FALSE)	
			}
			sortValues <- sort(values, na.last = NA)
			if (!identical(sortValues, values)) {
				stop(paste("In model", 
					omxQuotes(modelName),
					"the thresholds in column",
					omxQuotes(tName), "of matrix/algebra",
					omxQuotes(simplifyName(threshName, modelName)),
					"are not in ascending order.",
					"The current order is: ",
					omxQuotes(values), "and",
					"ascending order is: ",
					omxQuotes(sortValues),".",
					"Only the first", expectedThreshCount,
					"element(s) of this column are inspected."), call. = FALSE)	
			}
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
	if (is.data.frame(x)) {
		if (is.list(levels)) {
			return(data.frame(mapply(factor, x, levels, labels,
				MoreArgs=list(exclude = exclude, ordered = ordered), SIMPLIFY=FALSE),
				check.names = FALSE))
		} else {
			return(data.frame(lapply(x, factor, levels, labels, exclude, ordered),
				check.names = FALSE)) 
		}
	} else if (is.matrix(x)) {
		stop(paste("argument 'x' to mxFactor()",
		"is of illegal type matrix,",
		"legal types are vectors or data.frames"))
	} else {
		return(factor(x, levels, labels, exclude, ordered))
	}
}
