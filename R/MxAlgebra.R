#
#   Copyright 2007-2012 The OpenMx Project
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
setClassUnion("MxAlgebraFormula", c("call", "name", "language", "logical", "numeric"))

setClass(Class = "MxAlgebra",
	representation = representation(
		formula = "MxAlgebraFormula",
		name = "character",
		dirty = "logical",
		.dimnames = "MxListOrNull",
		initial = "matrix",
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



mxAlgebra <- function(expression, name = NA, dimnames = NA) {
	if (single.na(name)) {
		name <- imxUntitledName()
	}
	imxVerifyName(name, 0)
	retval <- new("MxAlgebra", NA, name)
	retval@formula <- match.call()$expression
    algebraErrorChecking(retval@formula, "mxAlgebra")
	if(!(length(dimnames) == 1 && is.na(dimnames))) {
		dimnames(retval) <- dimnames
	}
	return(retval)	
}

algebraErrorChecking <- function(formula, context) {
	if(length(formula) < 2) {
		return()
	}
	operator <- as.character(formula[[1]])
	if (identical(operator, "(")) {
	} else if (operator %in% omxSymbolTable[,'R.name']) {
	} else {
		msg <- paste("Unknown matrix operator or function",
			omxQuotes(operator), "in",
			deparse(width.cutoff = 400L, imxLocateFunction(context)))
		stop(msg, call. = FALSE)
	}
	lapply(formula[-1], algebraErrorChecking)
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

algebraSymbolCheck <- function(formula, name) {
	formula <- unlist(formulaList(formula))
	test <- sapply(formula, function(x) {!is.numeric(x) && !is.character(x)})
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
	matrixNumbers <- as.list(as.integer(-1 : (-length(matrixNames))))
	algebraNumbers <- as.list(as.integer(0 : (length(algebraNames) - 1)))
	names(matrixNumbers) <- matrixNames
	names(algebraNumbers) <- algebraNames
	retval <- eval(substitute(substitute(e, matrixNumbers), list(e = retval)))
	retval <- eval(substitute(substitute(e, algebraNumbers), list(e = retval)))
	retval <- substituteOperators(as.list(retval), algebra@name)
	algebraSymbolCheck(retval, algebra@name)
	return(list(algebra@initial, retval))
}

substituteOperators <- function(algebra, name) {
	if ((length(algebra) == 1) && (is.list(algebra))) {
		algebra <- list(0L, algebra[[1]])
	} else if ((length(algebra) > 1) && (!is.numeric(algebra[[1]]))) {
		if (as.character(algebra[[1]]) == '(') {
			return(substituteOperators(as.list(algebra[[2]]), name))
		}
		names <- omxSymbolTable["R.name"] == as.character(algebra[[1]])
        variableSymbols <- omxSymbolTable["Number.of.arguments"] == -1
		result <- omxSymbolTable[names & variableSymbols, "Num"]
		if (length(result) > 1) {
				msg <- paste("In algebra", omxQuotes(name), 
					"ambiguous function with name", algebra[[1]],
					"and", (length(algebra) - 1), "arguments.")
				stop(msg, call. = FALSE)
		} else if(length(result) == 1) {
			head <- as.integer(result[[1]])
			tail <- lapply(algebra[-1], substituteOperators, name)
			result <- append(tail, head, after=0)
			return(result)
		} else {
			length <- omxSymbolTable["Number.of.arguments"] == (length(algebra) - 1)
			result <- omxSymbolTable[names & length, "Num"]
			if (length(result) == 0) {
				msg <- paste("In algebra", omxQuotes(name),
					"could not find function with name", algebra[[1]],
					"and", (length(algebra) - 1), "arguments.")
				stop(msg, call. = FALSE)
			} else if (length(result) > 1) {
				msg <- paste("In algebra", omxQuotes(name),
					"ambiguous function with name", algebra[[1]],
					"and", (length(algebra) - 1), "arguments.")
				stop(msg, call. = FALSE)
			} else {
				head <- as.integer(result[[1]])
			    tail <- lapply(algebra[-1], substituteOperators, name)
				result <- append(tail, head, after=0)
				return(result)
			}
    	}
	}
	return(algebra)
}

checkEvaluation <- function(model, flatModel) {
	labelsData <- imxGenerateLabels(model)
	cache <- list()
	cache <- checkMatrixEvaluation(model, flatModel, labelsData, cache)
	tuple <- checkAlgebraEvaluation(model, flatModel, labelsData, cache)
	flatModel <- tuple[[1]]
	cache <- tuple[[2]]
	checkConstraintEvaluation(model, flatModel, labelsData, cache)
	return(flatModel)
}

checkMatrixEvaluation <- function(model, flatModel, labelsData, cache) {
	if(length(flatModel@matrices) == 0) { return(cache) }
	for(i in 1:length(flatModel@matrices)) {
		matrix <- flatModel@matrices[[i]]
		tuple <- evaluateMxObject(matrix@name, flatModel, labelsData, cache)
		cache <- tuple[[2]]
	}
	return(cache)
}

checkAlgebraEvaluation <- function(model, flatModel, labelsData, cache) {
	if(length(flatModel@algebras) == 0) { return(list(flatModel, cache)) }
	for(i in 1:length(flatModel@algebras)) {
		algebra <- flatModel@algebras[[i]]
		tuple <- evaluateMxObject(algebra@name, flatModel, labelsData, cache)
		result <- tuple[[1]]
		cache <- tuple[[2]]
		algebra <- flatModel[[algebra@name]]
		algebra@initial <- as.matrix(result)
		flatModel[[algebra@name]] <- algebra
	}
	return(list(flatModel, cache))
}

checkConstraintEvaluation <- function(model, flatModel, labelsData, cache) {
	if(length(flatModel@constraints) == 0) { return() }
	for(i in 1:length(flatModel@constraints)) {
		constraint <- flatModel@constraints[[i]]		
		lhsContext <- paste("the left-hand side of constraint", omxQuotes(constraint@name))
		rhsContext <- paste("the right-hand side of constraint", omxQuotes(constraint@name))		
		tuple <- evaluateAlgebraWithContext(flatModel[[constraint@alg1]], lhsContext, flatModel, labelsData, cache)
		lhs <- tuple[[1]]
		cache <- tuple[[2]]
		tuple <- evaluateAlgebraWithContext(flatModel[[constraint@alg2]], rhsContext, flatModel, labelsData, cache)
		rhs <- tuple[[1]]
		cache <- tuple[[2]]
		if (!all(dim(lhs) == dim(rhs))) {
			lhsName <- constraint@formula[[2]]
			rhsName <- constraint@formula[[3]]
			if (length(lhsName) == 1) {
				lhsName <- simplifyName(deparse(lhsName), model@name)
			} else {
				lhsName <- deparse(lhsName)
			}
			if (length(rhsName) == 1) {
				rhsName <- simplifyName(deparse(rhsName), model@name)
			} else {
				rhsName <- deparse(rhsName)
			}
			stop(paste("The algebras/matrices", 
				omxQuotes(c(lhsName, rhsName)),
				"in model", omxQuotes(model@name),
				"are in constraint", omxQuotes(simplifyName(constraint@name, model@name)),
				"and are not of identical dimensions. The left-hand side is",
				nrow(lhs), "x", ncol(lhs), "and the right-hand side is",
				nrow(rhs), "x", paste(ncol(rhs), ".", sep = '')), call. = FALSE)
		}
	}		
}

displayAlgebra <- function(mxAlgebra) {
	cat("mxAlgebra", omxQuotes(mxAlgebra@name), '\n')
	cat("@formula: ", deparse(mxAlgebra@formula, width.cutoff=500L), '\n')
	if (length(mxAlgebra@result) == 0) {
		cat("@result: (not yet computed) ")
	} else {
		cat("@result:\n")
	}
	print(mxAlgebra@result)
	if (is.null(dimnames(mxAlgebra))) {
			cat("dimnames: NULL\n")
	} else {
		cat("dimnames:\n")
		print(dimnames(mxAlgebra))
	}
}

setMethod("print", "MxAlgebra", function(x,...) { displayAlgebra(x) })
setMethod("show", "MxAlgebra", function(object) { displayAlgebra(object) })
