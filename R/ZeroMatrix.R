setClass(Class = "ZeroMatrix",
	representation = representation(),
	contains = "MxNonSymmetricMatrix")

setMethod("initialize", "ZeroMatrix",
	function(.Object, name, values, specification, nrow, ncol, byrow, free) {
		if (!single.na(values)) {
			warning("Ignoring values matrix for ZeroMatrix construction")
		}
		if (!single.na(specification)) {
			warning("Ignoring specification matrix for ZeroMatrix construction")
		}
		if (free) {
			warning("Ignoring \'free\' parameter for ZeroMatrix construction")
		}
		specification <- new("MxSparseMatrix", 0, nrow, ncol)
		values <- Matrix(0, nrow, ncol)
		return(callNextMethod(.Object, specification, values, name))
	}
)

setMethod("omxVerifyMatrix", "ZeroMatrix",
	function(.Object) {
		callNextMethod(.Object)	
		if(nnzero(.Object@specification) > 0) { 
			stop(paste("Specification matrix of zero matrix",
				omxQuotes(.Object@name), "is not empty."), call. = FALSE) 
		} 
		if(nnzero(.Object@values) > 0) { 
			stop(paste("Values matrix of zero matrix",
				omxQuotes(.Object@name), "is not empty"), call. = FALSE)
		} 
	}
)
