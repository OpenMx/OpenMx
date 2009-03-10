setClass(Class = "IdenMatrix",
	representation = representation(),
	contains = "MxNonSymmetricMatrix")
	
setMethod("initialize", "IdenMatrix",
	function(.Object, name, values, specification, nrow, ncol, byrow, free) {
		if (!single.na(values)) {
			warning("Ignoring values matrix for IdenMatrix construction")
		}
		if (!single.na(specification)) {
			warning("Ignoring specification matrix for IdenMatrix construction")
		}
		if (free) {
			warning("Ignoring \'free\' parameter for IdenMatrix construction")
		}		
		if (nrow != ncol) {
			stop("Non-square matrix attempted for IdenMatrix constructor")
		}
		specification <- new("MxSparseMatrix", 0, nrow, ncol)
		values <- Matrix(diag(nrow))
		return(callNextMethod(.Object, specification, values, name))
	}
)

setMethod("omxVerifyMatrix", "IdenMatrix",
	function(.Object) {
		callNextMethod(.Object)
		verifySquare(.Object)
		if(nnzero(.Object@specification) > 0) { stop("Specification matrix is not empty") } 
		if(!suppressWarnings(all(.Object@values == diag(nrow(.Object@values))))) 
			{ stop("Values matrix is not the identity matrix") }
	}
)
