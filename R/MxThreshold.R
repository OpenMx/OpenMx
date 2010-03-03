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

matchSingleDefinitionName <- function(targetName, threshNames, flatModel) {
	for(i in 1:length(threshNames)) {
		entity <- flatModel[[threshNames[[i]]]]
		entityCols <- dimnames(entity)[[2]]
		if (targetName %in% entityCols) {
			return(threshNames[[i]])
		}
	}
	return(as.integer(NA))
}

generateThresholdColumns <- function(flatModel, threshNames, dataName) {
	retval <- list()
	definitionNames <- dimnames(flatModel@datasets[[dataName]]@observed)[[2]]
	if (single.na(threshNames) || (length(threshNames) == 0)) {
		return(as.list(replicate(length(definitionNames), as.integer(NA))))
	}
	for(i in 1:length(definitionNames)) {
		targetName <- definitionNames[[i]]
		threshMatch <- matchSingleDefinitionName(targetName, threshNames, flatModel)
		if (is.na(threshMatch)) {
			retval[[i]] <- threshMatch
		} else {
			matchIndex <- omxLocateIndex(flatModel, threshMatch, flatModel@name)
			entity <- flatModel[[threshMatch]]
			entityCols <- dimnames(entity)[[2]]
			colIndex <- match(targetName, entityCols)
			retval[[i]] <- c(matchIndex, colIndex - 1L)
		}
	}
	return(retval)
}

convertThresholds <- function(flatModel, model, dataName, threshNames) {
	if (single.na(threshNames) || length(threshNames) == 0) {
		observed <- flatModel@datasets[[dataName]]@observed
		return(observed)
	}
	for(i in 1:length(threshNames)) {
		threshName <- threshNames[[i]]
		thresholds <- flatModel[[threshName]]
		if (is.null(thresholds)) {
			stop(paste("The thresholds matrix/algebra", omxQuotes(threshName), 
				"for model", omxQuotes(flatModel@name), 
				"does not exist."), call. = FALSE)
		} else if (is(thresholds, "MxAlgebra") || is(thresholds, "MxMatrix")) {
			flatModel <- convertSingleThreshold(flatModel, dataName, threshName)
		} else {
			stop(paste("The thresholds entity", omxQuotes(threshName), 
				"for model", omxQuotes(flatModel@name), 
				"is not a MxMatrix or MxAlgebra."), call. = FALSE)
		}
	}
	if (length(threshNames) == 1) {
		checkSingleThresholdEntity(flatModel, model, dataName, threshNames[[1]])
	}
	observed <- flatModel@datasets[[dataName]]@observed
	return(observed)
}

checkSingleThresholdEntity <- function(flatModel, model, dataName, threshName) {
	observed <- flatModel@datasets[[dataName]]@observed
	modelName <- flatModel@name
	thresholds <- eval(substitute(mxEval(x, model, compute=TRUE),
		list(x = as.symbol(threshName))))
	threshNames <- dimnames(thresholds)[[2]]
	for(i in 1:length(threshNames)) {
		tName <- threshNames[[i]]
		column <- thresholds[,i]
		count <- sum(!is.na(column))
		if (count != (length(levels(observed[,tName])) - 1)) {
			stop(paste("The number of thresholds in column",
				omxQuotes(threshNames[[i]]),
				"is not one less than the number of levels",
				"in model", 
				omxQuotes(modelName)), call. = FALSE)
		}
		values <- column[1:count]
		if (any(is.na(values))) {
			stop(paste("The thresholds in column",
				omxQuotes(tName),
				"contain NA values in between non-NA values",
				"in model",
				omxQuotes(modelName)), call. = FALSE)			
		}
		if (count < length(column) && 
				any(!is.na(column[count + 1:length(column)]))) {
			stop(paste("The thresholds in column",
				omxQuotes(tName),
				"contain NA values in between non-NA values",
				"in model",
				omxQuotes(modelName)), call. = FALSE)
		}
		sortValues <- sort(values)
		if (!all(sortValues == values)) {
			stop(paste("The thresholds in column",
				omxQuotes(tName),
				"are not in sorted order",
				"in model",
				omxQuotes(modelName)), call. = FALSE)	
		}
	}
}

convertSingleThreshold <- function(flatModel, dataName, threshName) {
	observed <- flatModel@datasets[[dataName]]@observed
	thresholds <- flatModel[[threshName]]
	modelName <- flatModel@name
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
		observed[,tName] <- as.ordered(observed[,tName])
	}
	flatModel@datasets[[dataName]]@observed <- observed
	return(flatModel)
}