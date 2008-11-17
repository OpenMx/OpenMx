setClass(Class = "FullMatrix",
	representation = representation(),
	contains = "MxNonSymmetricMatrix")

setMethod("initialize", "FullMatrix",
	function(.Object, name, values, specification, nrow, ncol, byrow, free) {
	    if (single.na(specification) && free) {
			specification <- new("MxSparseMatrix", matrix(NA, nrow, ncol))
	    } else if (single.na(specification)) {
			specification <- new("MxSparseMatrix", 0, nrow, ncol)
	    } else {
	    	specification <- new("MxSparseMatrix", specification, nrow, ncol)
	    }
	    if (single.na(values)) {
	    	values <- Matrix(0, nrow, ncol)
	    } else {
	    	values <- Matrix(values, nrow, ncol)
	    } 
		retval <- callNextMethod(.Object, specification, values, name) 
		return(retval)
	}
)

setMethod("verify", "FullMatrix",
	function(.Object) {
		callNextMethod(.Object)
	}
)
