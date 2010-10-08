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

definitionStartingValue <- function(defName, matrixName, flatModel, defvar.row = 1) {
	components <- unlist(strsplit(defName, omxSeparatorChar, fixed = TRUE))
	dataName <- paste(components[[1]], components[[2]], sep = '.')
	dataSet <- flatModel@datasets[[dataName]]
	if(is.null(dataSet)) {
		stop(paste("Could not find the dataset", 
			omxQuotes(simplifyName(dataName, flatModel@name)),
			"in model", omxQuotes(flatModel@name), 
			"used by the definition variable", 
			omxQuotes(simplifyName(defName, flatModel@name))),
			call. = FALSE)
	}
	if (!(components[[3]] %in% colnames(dataSet@observed))) {
		stop(paste("Could not find the column name", omxQuotes(components[[3]]),
			"in the dataset", omxQuotes(simplifyName(dataName, flatModel@name)), 
			"used by the definition variable", 
			omxQuotes(simplifyName(defName, flatModel@name))),
			call. = FALSE)
	}
	if (defvar.row < 1 || defvar.row > nrow(dataSet@observed)) {
		stop(paste("Row number", omxQuotes(defvar.row),
			"is out of bounds for definition variable",
			omxQuotes(simplifyName(defName, flatModel@name)), 
			"used in the context of",
			omxQuotes(simplifyName(matrixName, flatModel@name))), call. = FALSE)
	}
	return(dataSet@observed[defvar.row, components[[3]]])
}

populateDefInitialValues <- function(flatModel) {
	matrices <- flatModel@matrices
	if (length(matrices) > 0) {
		for(i in 1:length(matrices)) {
			flatModel@matrices[[i]]@values <- populateDefVarMatrix(matrices[[i]], flatModel)
		}
	}
	return(flatModel)
}

populateDefVarMatrix <- function(matrix, model, defvar.row = 1) {
	labels <- matrix@labels
	select <- !apply(labels, c(1,2), is.na) & apply(labels, c(1,2), omxIsDefinitionVariable)
	if (all(!select)) { return(matrix@values) }
	rows <- row(labels)[select]
	cols <- col(labels)[select]
	subs <- labels[select]
	value <- matrix@values
	for(i in 1:length(subs)) {
		startValue <- definitionStartingValue(subs[[i]], matrix@name, model, defvar.row)
		value[rows[[i]],cols[[i]]] <- startValue
	}
	return(value)
}

defVariableIsMatch <- function(defName, dataName) {
	components <- unlist(strsplit(defName, omxSeparatorChar, fixed = TRUE))
	target <- paste(components[[1]], components[[2]], sep = '.')
	return(identical(target, dataName))
}

omxFilterDefinitionVariables <- function(defVars, dataName) {
	if (length(defVars) == 0) return(defVars)
	dvNames <- names(defVars)	
	filter <- sapply(dvNames, defVariableIsMatch, dataName)
	return(defVars[filter])
}



