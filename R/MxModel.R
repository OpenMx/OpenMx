setClass(Class = "MxModel",
	representation = representation(
		matrices = "list",
		algebras = "list"))
		
setMethod("initialize", "MxModel",
	function(.Object, matrices=list(), algebras=list()) {
		.Object@matrices = matrices
		.Object@algebras = algebras
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

omxGenerateParameterList <- function(mxModel) {
	result <- list()
	for(i in 1:length(mxModel@matrices)) {
		result <- generaterParameterListHelper(mxModel@matrices[[i]], result, i - 1)
	}	
	return(result)
}
