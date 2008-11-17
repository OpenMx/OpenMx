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

setMethod("verify", "ZeroMatrix",
	function(.Object) {
		callNextMethod(.Object)	
		if(nnzero(.Object@specification) > 0) { stop("Specification matrix is not empty") } 
		if(nnzero(.Object@values) > 0) { stop("Values matrix is not empty") } 
	}
)
