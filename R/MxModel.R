setClass(Class = "MxModel",
	representation = representation(
		matrices = "list",
		algebras = "list",
		data = "data.frame"))
		
setMethod("initialize", "MxModel",
	function(.Object, matrices = list(), algebras = list(), data = data.frame()) {
		.Object@matrices = matrices
		.Object@algebras = algebras
		.Object@data = data
		return(.Object)
	}
)

setMethod("[[", "MxModel",
	function(x, i, j, ..., drop = FALSE) {
		first <- x@matrices[[i]]
		second <- x@algebras[[i]]
		if (is.null(first)) {
			return(second)
		} else {
			return(first)
		}	
	}
)

setReplaceMethod("[[", "MxModel",
	function(x, i, j, value) {
		if (is(value,"MxMatrix")) {
			if (!is.null(x@algebras[[i]])) {
				stop(paste(i, "is already an MxAlgebra object"))
			}
			x@matrices[[i]] <- value
		} else if (is(value,"MxAlgebra")) {
			if (!is.null(x@matrices[[i]])) {
				stop(paste(i, "is already an MxMatrix object"))
			}
			x@algebras[[i]] <- value		
		} else {
			stop(paste("Unknown type of value", value))
		}
		return(x)
	}
)

omxGenerateMatrixList <- function(mxModel) {
	return(lapply(mxModel@matrices, generateMatrixListHelper))
}

omxGenerateSimpleMatrixList <- function(mxModel) {
	retval <- lapply(mxModel@matrices, generateMatrixListHelper)
	return(lapply(retval, as.matrix))
}

omxGenerateAlgebraList <- function(mxModel) {
    mList <- omxGenerateMatrixList(mxModel)
    retval <- lapply(mxModel@algebras, generateAlgebraHelper, names(mList))
    return(retval)
}

omxGenerateParameterList <- function(mxModel) {
	result <- list()
	for(i in 1:length(mxModel@matrices)) {
		result <- generaterParameterListHelper(mxModel@matrices[[i]], result, i - 1)
	}	
	return(result)
}

omxGenerateValueList <- function(mxModel) {
	mList <- omxGenerateMatrixList(mxModel)
	pList <- omxGenerateParameterList(mxModel)
	retval <- vector()
	for(i in 1:length(pList)) {
		parameter <- pList[[i]]
		if (length(parameter) > 1) {
			values <- sapply(parameter, generateValueHelper, mList)
			if (!all(values == values[[1]])) {
				warning(paste('Parameter',i,'has multiple start values.',
					'Selecting', values[[1]]))
			}
			retval[i] <- values[[1]]
		} else {
			retval[i] <- generateValueHelper(parameter[[1]], mList)
		}
    }
	return(retval)	
}

options(stringsAsFactors = FALSE)
data(omxSymbolTable)
options(stringsAsFactors = TRUE)

generateAlgebraHelper <- function(algebra, lnames) {
	retval <- algebra@formula
	numbers <- as.list(as.double(-1 : (-length(lnames))))
	names(numbers) <- lnames
	retval <- eval(substitute(substitute(e, numbers), list(e = retval)))
	retval <- substituteOperators(retval)
	return(retval)
}

substituteOperators <- function(algebra) {
	if ((length(algebra) > 1) && (!is.numeric(algebra[[1]]))) {
		names <- omxSymbolTable["R.name"] == as.character(algebra[[1]])
        variableSymbols <- omxSymbolTable["Number.of.arguments"] == -1
		result <- omxSymbolTable[names & variableSymbols, "Num"]
		if (length(result) > 1) {
				stop(paste("Ambiguous function with name", algebra[[1]],
					"and", (length(algebra) - 1), "arguments"))
		} else if(length(result) == 1) {
			head <- as.double(result[[1]])
			tail <- lapply(algebra[-1], substituteOperators)
			return(list(head, tail))
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

generateValueHelper <- function(triple, mList) {
	mat <- triple[1] + 1
	row <- triple[2]
	col <- triple[3]
	return(mList[[mat]][row,col])
}

omxUpdateModelValues <- function(mxModel, values) {
	pList <- omxGenerateParameterList(mxModel)
	if(length(pList) != length(values)) {
		stop(paste("This model has", length(pList), 
			"parameters, but you have given me", length(values),
			"values"))
	}
	for(i in 1:length(pList)) {
		mxModel <- updateModelValueHelper(pList[[i]], values[[i]], mxModel)
    }
	return(mxModel)
}

updateModelValueHelper <- function(triples, value, mxModel) {
	for(i in 1:length(triples)) {
		triple <- triples[[i]]
		mat <- triple[1] + 1
		row <- triple[2]
		col <- triple[3]
		mxModel@matrices[[mat]]@values[row,col] <- value			
	}
	return(mxModel)
}
