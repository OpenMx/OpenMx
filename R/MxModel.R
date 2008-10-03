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