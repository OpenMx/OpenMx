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


omxQuotes <- function(name) {
	listTerms <- sapply(name, function(x) {paste("'", x, "'", sep = '')} )
	return(paste(listTerms, collapse=', '))
}

printOptions <- function(options) {
	retval <- ""
	select <- list()
	if (length(options) > 0) {
		for(i in 1:length(options)) {
			key <- names(options)[[i]] 
			if(!(key %in% names(getOption('mxOptimizerOptions')))) {
				select[[key]] <- options[[key]]
			} else if (options[[key]] != getOption('mxOptimizerOptions')[[key]]) {
				select[[key]] <- options[[key]]
			}
		}
	}	
	
	missingKeys <- setdiff(names(getOption('mxOptimizerOptions')), names(options))
	if (length(missingKeys) > 0) {
		for(i in 1:length(missingKeys)) {
			key <- missingKeys[[i]]
			select[[key]] <- "NULL"
		}
	}
	
	if (length(select) == 0) {
		return(retval)
	}
	for(i in 1:length(select)) {
		key <- names(select)[[i]]
		value <- select[[i]]
		retval <- paste(retval, omxQuotes(key), '=',
				omxQuotes(value))
		if (i < length(select)) {
			retval <- paste(retval, ',', sep='')
		}
	}
	return(retval)
}

displayModel <- function(model, expand = FALSE) {
	cat("MxModel", omxQuotes(model@name), '\n')
	cat("type :", omxTypeName(model), '\n')
	cat("matrices :", omxQuotes(names(model@matrices)), '\n')
	cat("algebras :", omxQuotes(names(model@algebras)), '\n')
	cat("constraints :", omxQuotes(names(model@constraints)), '\n')
	if (length(model@latentVars) == 0) {
		cat("latentVars : none\n")
	} else {
		cat("latentVars :", omxQuotes(model@latentVars), '\n')
	}
	if (length(model@manifestVars) == 0) {
		cat("manifestVars : none\n")
	} else {
		cat("manifestVars :", omxQuotes(model@manifestVars), '\n')
	}
	data <- model@data
	if (is.null(data)) {
		cat("data : NULL\n")
	} else {
		cat("data matrix :", nrow(data@matrix), 
			"x", ncol(data@matrix), "\n")
		if(length(data@means) == 1 && is.na(data@means)) {
			cat("data means : NA\n")
		} else {
			cat("data means : 1 x", length(data@means), "\n")
		}
		cat("data type:", omxQuotes(data@type), '\n')
	}
	cat("submodels :", omxQuotes(names(model@submodels)), '\n')
	objective <- model@objective
	if (is.null(objective)) {
		objectiveType <- "NULL"
	} else {
		objectiveType <- class(objective)[[1]]
	}
	cat("objective :", objectiveType, '\n')
	cat("independent :", model@independent, '\n')
	cat("options :", printOptions(model@options), '\n')
	cat("output :", length(model@output) > 0, '\n')
	if(expand) {
		if(length(model@matrices) > 0) {
			cat("\n--------MATRICES--------\n")
			lapply(model@matrices, print)
		}
		if(length(model@algebras) > 0) {
			cat("\n--------ALGEBRAS--------\n")
			lapply(model@algebras, print)
		}
		if(length(model@constraints) > 0) {
			cat("\n--------CONSTRAINTS--------\n")
			lapply(model@constraints, print)
		}
		if(!is.null(model@data) > 0) {
			cat("\n--------DATA--------\n")
			print(model@data)
		}
		if(!is.null(model@objective) > 0) {
			cat("\n--------OBJECTIVE FUNCTION--------\n")
			print(model@objective)
		}		
		if(length(model@output) > 0) {
			cat("\n--------OUTPUT--------\n")
			print(model@output)
		}
		if(length(model@submodels) > 0) {
			cat("\n--------SUBMODELS--------\n")
			lapply(model@submodels, print)
		}
		if(length(model@options) > 0) {
			cat("\n--------OPTIONS--------\n")
			print(model@options)
		}
	}
	invisible(model)
}

setMethod("print", "MxModel", function(x,...) { 
	displayModel(x, TRUE) 
})

setMethod("show", "MxModel", function(object) { 
	displayModel(object) 
})

