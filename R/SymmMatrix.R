setClass(Class = "SymmMatrix",
	representation = representation(),
	contains = "MxSymmetricMatrix")
	
setMethod("initialize", "SymmMatrix",
	function(.Object, name, values, specification, nrow, ncol, byrow, free) {
		if (nrow != ncol) {
			stop("Non-square matrix attempted for SymmMatrix constructor")
		}
		if (!single.na(values) && is.vector(values)) {
		    mvalues <- matrix(0, nrow, ncol)
		    mvalues[lower.tri(mvalues, diag = TRUE)] <- values
		    mvalues <- mvalues + t(mvalues) - diag(mvalues) * diag(nrow)
		    values <- mvalues
		}
		if (!single.na(specification) && is.vector(specification)) {
		    mspec <- matrix(0, nrow, ncol)
		    mspec[lower.tri(mspec, diag = TRUE)] <- specification
		    mspec <- mspec + t(mspec) - diag(mspec) * diag(nrow)
		    specification <- mspec
		}				
	    if (single.na(specification) && free) {
			specification <- new("MxSymmetricSparse", matrix(NA, nrow, ncol))
	    } else if (single.na(specification)){
			specification <- new("MxSymmetricSparse", 0, nrow, ncol)
	    } else {
	    	specification <- new("MxSymmetricSparse", matrix(specification, nrow, ncol))
	    }	    
	    if (single.na(values)) {
	    	values <- new("MxSymmetricSparse", matrix(0, nrow, ncol))
	    } else {
	    	values <- new("MxSymmetricSparse", matrix(values, nrow, ncol))
	    }
		retval <- callNextMethod(.Object, specification, values, name)
		return(retval)
	}
)

setMethod("verify", "SymmMatrix",
	function(.Object) {
		callNextMethod(.Object)
		verifySquare(.Object)
		values <- .Object@values
		if (!all(values == t(values))) {
			stop(paste("Symmetric matrix",.Object@name,"is not symmetric!"))
		}
	}
)
