setClassUnion("MxAlgebraFormula", c("call", "name", "logical"))

setClass(Class = "MxAlgebra",
	representation = representation(
		formula = "MxAlgebraFormula",
		name = "character",
		result = "Matrix"))
		
setMethod("initialize", "MxAlgebra",
	function(.Object, formula = NA, name = NA) {
		.Object@formula <- sys.call(which=-3)[[3]]
		if (is.na(name)) {
			.Object@name <- omxUntitledName()
		} else {
			.Object@name <- name
		}
		return(.Object)
	}
)

mxAlgebra <- function(expression, name = NA) {
	retval <- new("MxAlgebra", NA, name)
	retval@formula <- match.call()$expression
	return(retval)	
}

options(stringsAsFactors = FALSE)
data(omxSymbolTable)
options(stringsAsFactors = TRUE)

generateAlgebraHelper <- function(algebra, matrixNames, algebraNames) {
	retval <- algebra@formula
	matrixNumbers <- as.list(as.double(-1 : (-length(matrixNames))))
	algebraNumbers <- as.list(as.double(0 : (length(algebraNames) - 1)))
	names(matrixNumbers) <- matrixNames
	names(algebraNumbers) <- algebraNames
	retval <- eval(substitute(substitute(e, matrixNumbers), list(e = retval)))
	retval <- eval(substitute(substitute(e, algebraNumbers), list(e = retval)))
	retval <- substituteOperators(as.list(retval))
	return(retval)
}

substituteOperators <- function(algebra) {
	if ((length(algebra) == 1) && (is.list(algebra))) {
		algebra <- list(0, algebra[[1]])
	} else if ((length(algebra) > 1) && (!is.numeric(algebra[[1]]))) {
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
