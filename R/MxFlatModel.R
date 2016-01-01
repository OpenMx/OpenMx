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


##' MxFlatModel
##'
##' This is an internal class and should not be used.
##'
##' @aliases
##' MxFlatModel-class
##' $,MxFlatModel-method
##' $<-,MxFlatModel-method
##' [[,MxFlatModel-method
##' [[<-,MxFlatModel-method
##' print,MxFlatModel-method
##' show,MxFlatModel-method
##' names,MxFlatModel-method
setClass(Class = "MxFlatModel",
	representation = representation(
		expectations = "list",
		fitfunctions = "list",
		datasets = "list",
		constMatrices = "list",
		freeMatrices = "list",
		parameters = "list"
	),
	contains = "MxModel")
	
setMethod("initialize", "MxFlatModel",
	function(.Object, model) {
		.Object@name <- model@name
		return(.Object)
	}
)

setMethod("[[", "MxFlatModel",
	function(x, i, j, ..., drop = FALSE) {
		return(flatExtractMethod(x, i))
	}
)

setReplaceMethod("[[", "MxFlatModel",
	function(x, i, j, value) {
		return(flatReplaceMethod(x, i, value))
	}
)

setMethod("$", "MxFlatModel",
	function(x, name) {
		return(flatExtractMethod(x, name))
	}
)

setReplaceMethod("$", "MxFlatModel",
	function(x, name, value) {
		return(flatReplaceMethod(x, name, value))
	}
)

setMethod("names", "MxFlatModel", 
	function(x) {
			return(flatNamesMethod(x))
	}
)

flatNamesMethod <- function(model) {
	return(flatNamespaceList(model))
}

flatExtractMethod <- function(model, index) {
	return(flatNamespaceSearch(model, index))
}

flatReplaceMethod <- function(model, index, value) {
	return(flatNamespaceSearchReplace(model, index, value))
}

##' imxCheckVariables
##'
##' This is an internal function exported for those people who know
##' what they are doing.
##' 
##' @param flatModel flatModel
##' @param namespace namespace
imxCheckVariables <- function(flatModel, namespace) {
	datasets <- flatModel@datasets
	retval <- list()
	retval$parameters <- namespace$parameters
	retval$values <- namespace$values
	retval$startvals <- list()
	retval$bounds <- list()
	if(length(flatModel@matrices) > 0) {
		startVals <- list()
		freeVars <- list()
		fixedVars <- list()
		bounds <- list()
		varlocations <- list()
		for(i in 1:length(flatModel@matrices)) {
			result <- checkVariablesHelper(flatModel@matrices[[i]], startVals, 
				freeVars, fixedVars, bounds, varlocations, flatModel@name, datasets)
			startVals <- result[[1]]
			freeVars <- result[[2]]
			fixedVars <- result[[3]]
			bounds <- result[[4]]
			varlocations <- result[[5]]
		}
		retval$startvals <- startVals
		retval$bounds <- bounds
	}
	return(retval)
}

rowColToString <- function(row, col) {
	return(paste("(", row, ", ", col, ")", sep = ""))
}

checkVariablesHelper <- function(matrix, startVals, freeVars,
		fixedVars, bounds, varlocations, modelname, datasets) {
	labels <- matrix@labels
	free <- matrix@free
	values <- matrix@values
	lbounds <- matrix@lbound
	ubounds <- matrix@ubound
	select <- !is.na(labels)
	free <- free[select]
	values <- values[select]
	lbounds <- lbounds[select]
	ubounds <- ubounds[select]
	rows <- row(labels)[select]
	cols <- col(labels)[select]
	labels <- labels[select]
	if(length(labels) > 0) {
		for(i in 1:length(labels)) {
			label <- labels[[i]]
			isFree <- free[[i]]
			value <- values[[i]]
			lbound <- lbounds[[i]]
			ubound <- ubounds[[i]]
			row <- rows[[i]]
			col <- cols[[i]]
			if (imxIsDefinitionVariable(label)) {
				if (isFree) {
					stop(paste("The definition variable", omxQuotes(label),
						"has been assigned to a free parameter",
						"in matrix", omxQuotes(simplifyName(matrix@name, modelname))), call. = FALSE)
				}
				result <- strsplit(label, '.', fixed=TRUE)[[1]]
				if(length(result) != 3) {
					stop("Internal error: definition variable does not have three pieces", call. = FALSE)
				}
				dataname <- paste(result[[1]], '.', result[[2]], sep = '')
				targetdata <- datasets[[dataname]]
				if (is.null(targetdata)) {
					stop(paste("The definition variable", omxQuotes(label),
						"in matrix", omxQuotes(simplifyName(matrix@name, modelname)),
						"refers to a data set that does not exist"), call. = FALSE)
				}
				targetNames <- dimnames(targetdata@observed)
				if (is.null(targetNames) || is.null(targetNames[[2]])) {
					stop(paste("The definition variable", omxQuotes(label),
						"in matrix", omxQuotes(simplifyName(matrix@name, modelname)),
						"refers to a data set that does not contain column names"), call. = FALSE)	
				}
				if (!(result[[3]] %in% targetNames[[2]])) {
					stop(paste("The definition variable", omxQuotes(label),
						"in matrix", omxQuotes(simplifyName(matrix@name, modelname)),
						"refers to a data set that does not contain",
						"a column with name", omxQuotes(result[[3]])), call. = FALSE)	
				}
			} else if (isFree) {
				if (hasSquareBrackets(label)) {
					stop(paste("The label with square brackets",
						"has been assigned to a free parameter",
						"in matrix", omxQuotes(simplifyName(matrix@name, modelname)),
						"at row", row, "and column", col), call. = FALSE)
				} else if (label %in% fixedVars) {
					stop(paste("The label", omxQuotes(label),
						"has been assigned to a free parameter",
						"and a fixed value!"), call. = FALSE)
				} else if (label %in% freeVars && !is.na(value) && length(startVals[[label]]) &&
					   startVals[[label]] != value) {
					loc <- varlocations[[label]]
					stop(paste("The free parameter", omxQuotes(label),
						"has been assigned multiple starting values!",
						"See matrix",
						omxQuotes(simplifyName(matrix@name, modelname)), 
						"at location", 
						rowColToString(row, col),
						"and matrix", 
						omxQuotes(simplifyName(loc[[1]], modelname)), 
						"at location", 
						rowColToString(loc[[2]], loc[[3]]), 
						"If you want to randomly select one of these values, call",
						"model <- omxAssignFirstParameters(model) before running again."), call. = FALSE)
				} else if (label %in% freeVars && 
								!identicalNA(lbound, bounds[[label]][[1]])) {
					loc <- varlocations[[label]]
					stop(paste("The free parameter", omxQuotes(label),
						"has been assigned multiple lower bounds!",
						"See matrix",
						omxQuotes(simplifyName(matrix@name, modelname)), 
						"at location", 
						rowColToString(row, col),
						"and matrix", 
						omxQuotes(simplifyName(loc[[1]], modelname)), 
						"at location", 
						rowColToString(loc[[2]], loc[[3]])), call. = FALSE)
				} else if (label %in% freeVars && 
								!identicalNA(ubound, bounds[[label]][[2]])) {
					loc <- varlocations[[label]]
					stop(paste("The free parameter", omxQuotes(label),
						"has been assigned multiple upper bounds!",
						"See matrix",
						omxQuotes(simplifyName(matrix@name, modelname)), 
						"at location", 
						rowColToString(row, col),
						"and matrix", 
						omxQuotes(simplifyName(loc[[1]], modelname)), 
						"at location", 
						rowColToString(loc[[2]], loc[[3]])), call. = FALSE)
				} else {
					if (!is.na(value)) startVals[[label]] <- value
					freeVars <- union(freeVars, label)
					bounds[[label]] <- c(lbound, ubound)
				}
			} else {
				if (hasSquareBrackets(label)) {
				} else if (label %in% freeVars) {
					stop(paste("The label", omxQuotes(label),
						"has been assigned to a fixed value",
						"and a free parameter!"), call. = FALSE)
				} else if (label %in% fixedVars && !is.na(value) && length(startVals[[label]]) &&
					   startVals[[label]] != value) {
					loc <- varlocations[[label]]
					stop(paste("The fixed variable", omxQuotes(label),
						"has been assigned multiple starting values!",
						"See matrix",
						omxQuotes(simplifyName(matrix@name, modelname)), 
						"at location", 
						rowColToString(row, col),
						"and matrix", 
						omxQuotes(simplifyName(loc[[1]], modelname)), 
						"at location", 
						rowColToString(loc[[2]], loc[[3]]),
						"If you want to randomly select one of these values, call",
						"model <- omxAssignFirstParameters(model)",
						"before running again."), call. = FALSE)
				} else {
					if (!is.na(value)) startVals[[label]] <- value
					fixedVars <- union(fixedVars, label)
				}
			}
			if(is.null(varlocations[[label]])) {
				varlocations[[label]] <- c(matrix@name, row, col)
			}
		}
	}
	return(list(startVals, freeVars, fixedVars, bounds, varlocations))
}

identicalNA <- function(x, y) {
	return((is.na(x) && is.na(y)) || (identical(x,y)))
}

displayFlatModel <- function(fm) {
	cat("expectations : ")
	print(fm@expectations)
	cat("fitfunctions : ")
	print(fm@fitfunctions)
	cat("compute : ")
	print(fm@compute)
	cat("datasets :", length(fm@datasets), '\n') 
}

setMethod("print", "MxFlatModel", function(x,...) { callNextMethod(); displayFlatModel(x) })
setMethod("show", "MxFlatModel", function(object) { callNextMethod(); displayFlatModel(object) })
