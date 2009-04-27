#
#   Copyright 2007-2009 The OpenMx Project
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


setClass(Class = "MxFlatModel",
	representation = representation(
		objectives = "list",
		datasets = "list"
	),
	contains = "MxModel")
	
setMethod("initialize", "MxFlatModel",
	function(.Object, model, objectives, datasets) {
		modelSlotNames <- slotNames(model)
		for(i in 1:length(modelSlotNames)) {
			name <- modelSlotNames[[i]]
			slot(.Object, name) <- slot(model, name)
		}
		.Object@objectives <- objectives
		.Object@datasets <- datasets
		return(.Object)
	}
)	

omxGenerateDefinitionNames <- function(datasets) {
	nameList <- lapply(datasets, 
		function(x) { dimnames(x@matrix)[[2]] })
	result <- list()
	if(length(nameList) > 0) {
		for(i in 1:length(nameList)) {
			colNames <- nameList[[i]]
			if(length(colNames) > 0) {
				for(j in 1:length(colNames)) {
					name <- colNames[[j]]
					result[[name]] <- c(i - 1, j - 1)
				}
			}	
		}
	}
	return(result)
}

omxCheckFreeVariables <- function(model, defLocations) {
	if(length(model@matrices) > 0) {
		startVals <- list()
		freeVars <- list()
		fixedVars <- list()
		defNames <- names(defLocations)
		for(i in 1:length(model@matrices)) {
			result <- omxCheckFreeVariablesHelper(model@matrices[[i]], startVals, 
				freeVars, fixedVars, defNames)
			startVals <- result[[1]]
			freeVars <- result[[2]]
			fixedVars <- result[[3]]
		}
	}
}

omxCheckFreeVariablesHelper <- function(matrix, startVals, freeVars,
		fixedVars, defNames) {
	labels <- matrix@labels
	free <- matrix@free
	values <- matrix@values
	labels <- labels[!is.na(labels)]
	free <- free[!is.na(labels)]
	values <- values[!is.na(labels)]
	if(length(labels) > 0) {
		for(i in 1:length(labels)) {
			label <- labels[[i]]
			isFree <- free[[i]]
			value <- values[[i]]
			if (label %in% defNames) {
			} else if (isFree) {
				if (label %in% fixedVars) {
					stop(paste("The label", omxQuotes(label),
						"has been assigned to a free parameter",
						"and a fixed value!"), call. = FALSE)
				} else if (label %in% freeVars && startVals[[label]] != value) {
					stop(paste("The free parameter", omxQuotes(label),
						"has been assigned multiple starting values!"),
						 call. = FALSE)
				} else {
					startVals[[label]] <- value
					freeVars <- append(freeVars, label)
				}
			} else {
				if (label %in% freeVars) {
					stop(paste("The label", omxQuotes(label),
						"has been assigned to a free parameter",
						"and a fixed value!"), call. = FALSE)
				} else if (label %in% fixedVars && startVals[[label]] != value) {
					stop(paste("The fixed variable", omxQuotes(label),
						"has been assigned multiple starting values!"),
						 call. = FALSE)
				} else {
					startVals[[label]] <- value
					fixedVars <- append(fixedVars, label)
				}
			}
		}
	}
	return(list(startVals, freeVars, fixedVars))
}

setMethod("print", "MxFlatModel", function(x,...) {
	callNextMethod()
	cat("objectives : ")
	print(x@objectives)
	cat("datasets :", length(x@datasets), '\n') 
})

setMethod("show", "MxFlatModel", function(object) { 
	callNextMethod()
	cat("objectives : ")
	print(object@objectives)
	cat("datasets :", length(object@datasets), '\n') 
})
