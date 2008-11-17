setClass(Class = "UnitMatrix",
	representation = representation(),
	contains = "MxNonSymmetricMatrix")
	
setMethod("initialize", "UnitMatrix",
	function(.Object, name, values, specification, nrow, ncol, byrow, free) {
		if (!single.na(values)) {
			warning("Ignoring values matrix for UnitMatrix construction")
		}
		if (!single.na(specification)) {
			warning("Ignoring specification matrix for UnitMatrix construction")
		}
		if (free) {
			warning("Ignoring \'free\' parameter for UnitMatrix construction")
		}
		specification <- new("MxSparseMatrix", 0, nrow, ncol)
		values <- Matrix(1, nrow, ncol)
		return(callNextMethod(.Object, specification, values, name))
	}
)

setMethod("verify", "UnitMatrix",
	function(.Object) {
		callNextMethod(.Object)		
		if(nnzero(.Object@specification) > 0) { stop("Specification matrix is not empty") } 
		if(nnzero(.Object@values - 1) > 0) { stop("Values matrix has non unit entries") } 
	}
)
