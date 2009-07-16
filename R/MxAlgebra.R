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

setClassUnion("MxListOrNull", c("list", "NULL"))
setClassUnion("MxAlgebraFormula", c("call", "name", "logical"))

setClass(Class = "MxAlgebra",
	representation = representation(
		formula = "MxAlgebraFormula",
		name = "character",
		dirty = "logical",
		.dimnames = "MxListOrNull",
		result = "matrix"))
		
setMethod("initialize", "MxAlgebra",
	function(.Object, formula, name) {
		.Object@formula <- sys.call(which=-3)[[3]]
		.Object@name <- name
		.Object@dirty <- FALSE
		.Object@.dimnames <- NULL
		return(.Object)
	}
)

setMethod("dimnames", "MxAlgebra",
	function(x) { x@.dimnames }
)

setReplaceMethod("dimnames", "MxAlgebra",
	function(x, value) {
		if (is.null(value)) {
		} else if (!is.list(value)) {
			stop("dimnames of MxAlgebra object must be either NULL or list of length 2")
		} else if (length(value) != 2) {
			stop("dimnames of MxAlgebra object must be either NULL or list of length 2")
		}
		x@.dimnames <- value
		return(x)
	}
)

mxAlgebra <- function(expression, name = NA) {
	if (is.na(name)) {
		name <- omxUntitledName()
	}
	omxVerifyName(name)
	retval <- new("MxAlgebra", NA, name)
	retval@formula <- match.call()$expression
	return(retval)	
}

defStringsAsFactors <- getOption('stringsAsFactors')
options(stringsAsFactors = FALSE)
data(omxSymbolTable)
options(stringsAsFactors = defStringsAsFactors)


formulaList <- function(x) {
	retval <- as.list(x)
	retval <- lapply(retval, function(x) {
		if(is.call(x)) {formulaList(x)} else {x}
	})
	return(retval)
}

algebraNumericCheck <- function(formula, name) {
	formula <- unlist(formulaList(formula))
	test <- sapply(formula, is.numeric)
	if(any(test)) {
		msg <- paste("There is a numeric operand in",
			"the algebra named", omxQuotes(name))
		stop(msg, call. = FALSE)
	}
}

algebraSymbolCheck <- function(formula, name) {
	formula <- unlist(formulaList(formula))
	test <- sapply(formula, function(x) {!is.numeric(x)})
	if(length(formula[test]) == 1) {
		msg <- paste("The reference", omxQuotes(formula[test]),
			"is unknown in the algebra named", omxQuotes(name))
		stop(msg, call. = FALSE)
	} else if (length(formula[test]) > 1) {
		msg <- paste("The references", omxQuotes(formula[test]),
			"are unknown in the algebra named", omxQuotes(name))
		stop(msg, call. = FALSE)		
	}
}

generateAlgebraHelper <- function(algebra, matrixNames, algebraNames) {
	retval <- algebra@formula
	algebraNumericCheck(retval, algebra@name)
	matrixNumbers <- as.list(as.double(-1 : (-length(matrixNames))))
	algebraNumbers <- as.list(as.double(0 : (length(algebraNames) - 1)))
	names(matrixNumbers) <- matrixNames
	names(algebraNumbers) <- algebraNames
	retval <- eval(substitute(substitute(e, matrixNumbers), list(e = retval)))
	retval <- eval(substitute(substitute(e, algebraNumbers), list(e = retval)))
	retval <- substituteOperators(as.list(retval))
	algebraSymbolCheck(retval, algebra@name)
	return(retval)
}

substituteOperators <- function(algebra) {
	if ((length(algebra) == 1) && (is.list(algebra))) {
		algebra <- list(0, algebra[[1]])
	} else if ((length(algebra) > 1) && (!is.numeric(algebra[[1]]))) {
		if (as.character(algebra[[1]]) == '(') {
			return(substituteOperators(as.list(algebra[[2]])))
		}
		names <- omxSymbolTable["R.name"] == as.character(algebra[[1]])
        variableSymbols <- omxSymbolTable["Number.of.arguments"] == -1
		result <- omxSymbolTable[names & variableSymbols, "Num"]
		if (length(result) > 1) {
				stop(paste("Ambiguous function with name", algebra[[1]],
					"and", (length(algebra) - 1), "arguments"))
		} else if(length(result) == 1) {
			head <- as.double(result[[1]])
			tail <- lapply(algebra[-1], substituteOperators)
			result <- append(tail, head, after=0)
			return(result)
		} else {
			length <- omxSymbolTable["Number.of.arguments"] == (length(algebra) - 1)
			result <- omxSymbolTable[names & length, "Num"]
			if (length(result) == 0) {
				stop(paste("Could not find function with name", algebra[[1]],
					"and", (length(algebra) - 1), "arguments"))
			} else if (length(result) > 1) {
				stop(paste("Ambiguous function with name", algebra[[1]],
					"and", (length(algebra) - 1), "arguments"))
			} else {
				head <- as.double(result[[1]])
			    tail <- lapply(algebra[-1], substituteOperators)
				result <- append(tail, head, after=0)
				return(result)
			}
    	}
	}
	return(algebra)
}

displayAlgebra <- function(mxAlgebra) {
	cat("mxAlgebra", omxQuotes(mxAlgebra@name), '\n')
	cat("formula: ", deparse(mxAlgebra@formula, width.cutoff=500L), '\n')
	if (is.null(dimnames(mxAlgebra))) {
			cat("dimnames: NULL\n")
	} else {
		cat("dimnames:\n")
		print(dimnames(mxAlgebra))
	}
}

setMethod("print", "MxAlgebra", function(x,...) { displayAlgebra(x) })
setMethod("show", "MxAlgebra", function(object) { displayAlgebra(object) })