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

##' omxQuotes
##'
##' Quote helper function, often for error messages.
##' 
##' @param name a character vector
##' @details
##' This is a helper function for creating a nicely
##' put together formatted string.
##' @return
##' a character string
##' @examples
##' omxQuotes(c("Oh", "blah", "dee", "Oh", "blah", "da"))
##' omxQuotes(c("A", "S", "F"))
##' omxQuotes("Hello World")
omxQuotes <- function(name) {
	if (length(name) == 0) return('NULL')
	if (typeof(name) == "S4") return(omxQuotes(list(name)))
	listTerms <- sapply(name, function(x) {
		if (typeof(x) == "S4") {
			return(class(x))
		} else if (is.na(x)) {
			return(as.character(x))
		} else {
			paste("'", x, "'", sep = '')
		}
	} )
	if (length(listTerms) == 2) {
		return(paste(listTerms[1], ' and ', listTerms[2], sep = ''))
	} else if (length(listTerms) > 2) {
		return(paste(paste(listTerms[1:length(listTerms) - 1], collapse=', '),
			', and ', listTerms[[length(listTerms)]], sep = ''))
	} else {
		return(listTerms)
	}
}

printOptions <- function(options) {
	retval <- ""
	
	if (length(options) == 0) {
		return(retval)
	}
	for(i in 1:length(options)) {
		key <- names(options)[[i]]
		value <- options[[i]]
		retval <- paste(retval, omxQuotes(key), '=',
				omxQuotes(value))
		if (i < length(options)) {
			retval <- paste(retval, ',', sep='')
		}
	}
	return(retval)
}

displayModel <- function(model, expand = FALSE) {
	cat("MxModel", omxQuotes(model@name), '\n')
	cat("type :", imxTypeName(model), '\n')
	cat("$matrices :", omxQuotes(names(model@matrices)), '\n')
	cat("$algebras :", omxQuotes(names(model@algebras)), '\n')
	cat("$constraints :", omxQuotes(names(model@constraints)), '\n')
	cat("$intervals :", omxQuotes(names(model@intervals)), '\n')
	if (length(model@latentVars) == 0) {
		cat("$latentVars : none\n")
	} else if (is.character(model@latentVars)) {
		cat("$latentVars :", omxQuotes(model@latentVars), '\n')
	} else {
		cat("$latentVars :\n")
		print(format(model@latentVars))
	}
	if (length(model@manifestVars) == 0) {
		cat("$manifestVars : none\n")
	} else if (is.character(model@manifestVars)) {
		cat("$manifestVars :", omxQuotes(model@manifestVars), '\n')
	} else {
		cat("$manifestVars :\n")
		print(format(model@manifestVars))
	}
	data <- model@data
	if (is.null(data)) {
		cat("$data : NULL\n")
	} else {
		if (is(data, "MxDataDynamic")) {
			cat("$data type:", omxQuotes(data@type), '\n')
			cat("$data$expectation :", omxQuotes(data@expectation), "\n")
		} else {
			cat("$data :", nrow(data@observed), 
			    "x", ncol(data@observed), "\n")
			if(length(data@means) == 1 && is.na(data@means)) {
				cat("$data means : NA\n")
			} else {
				cat("$data means : 1 x", length(data@means), "\n")
			}
			cat("$data type:", omxQuotes(data@type), '\n')
		}
	}
	cat("$submodels :", omxQuotes(names(model@submodels)), '\n')
	expectation <- model@expectation
	fitfunction <- model@fitfunction
	compute <- model@compute
	if (is.null(expectation)) {
		expectationType <- "NULL"
	} else {
		expectationType <- class(expectation)[[1]]
	}
	cat("$expectation :", expectationType, '\n')

	if (is.null(fitfunction)) {
		fitfunctionType <- "NULL"
	} else {
		fitfunctionType <- class(fitfunction)[[1]]
	}
	cat("$fitfunction :", fitfunctionType, '\n')

	if (is.null(compute)) {
		computeType <- "NULL"
	} else {
		computeType <- class(compute)[[1]]
	}
	cat("$compute :", computeType, '\n')
	cat("$independent :", model@independent, '\n')
	cat("$options :", printOptions(model@options), '\n')
	cat("$output :", length(model@output) > 0, '\n')
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
		if(!is.null(model@expectation) > 0) {
			cat("\n--------EXPECTATION FUNCTION--------\n")
			print(model@expectation)
		}
		if(!is.null(model@fitfunction) > 0) {
			cat("\n--------FIT FUNCTION--------\n")
			print(model@fitfunction)
		}
		if(!is.null(model@compute) > 0) {
			cat("\n--------COMPUTE--------\n")
			print(model@compute)
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
