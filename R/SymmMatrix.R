setClass(Class = "SymmMatrix",
	representation = representation(),
	contains = "MxSymmetricMatrix")
	
setMethod("initialize", "SymmMatrix",
	function(.Object, name, values, specification, nrow, ncol, byrow, free) {
		if (nrow != ncol) {
			stop("Non-square matrix attempted for SymmMatrix constructor")
		}
		if (!single.na(values) && is.vector(values)) {
			if (length(values) == (nrow * (nrow + 1) / 2)) {
			    mvalues <- matrix(0, nrow, ncol)
			    mvalues[lower.tri(mvalues, diag = TRUE)] <- values
			    mvalues <- mvalues + t(mvalues) - diag(mvalues) * diag(nrow)
			    values <- mvalues
			} else if (length(values) == (nrow * ncol)) {
				values <- matrix(values, nrow, ncol)
			} else {
				stop("Invalid length of values matrix for SymmMatrix constructor")
			}
		}
		if (!single.na(specification) && is.vector(specification)) {
			if (length(specification) == (nrow * (nrow + 1) / 2)) {
			    mspec <- matrix(0, nrow, ncol)
			    mspec[lower.tri(mspec, diag = TRUE)] <- specification
			    mspec <- mspec + t(mspec) - diag(mspec) * diag(nrow)
			    specification <- mspec
			} else if (length(specification) == (nrow * ncol)) {
				specification <- matrix(specification, nrow, ncol)
			} else {
				stop("Invalid length of specification matrix for SymmMatrix constructor")
			}
		}
		if (is(specification, "MxSymmetricSparse")) {
		} else if (single.na(specification) && free) {
			specification <- new("MxSymmetricSparse", matrix(NA, nrow, ncol))
	    } else if (single.na(specification)){
			specification <- new("MxSymmetricSparse", 0, nrow, ncol)
		} else if (is(specification, "Matrix")) {
			specification <- new("MxSymmetricSparse", as.matrix(specification))
	    } else {
	    	specification <- new("MxSymmetricSparse", matrix(specification, nrow, ncol))
	    }
	    if (is(values, "MxSymmetricSparse")) {
	    } else if (single.na(values)) {
	    	values <- new("MxSymmetricSparse", matrix(0, nrow, ncol))
	    } else if (is(values, "Matrix")) {
			values <- new("MxSymmetricSparse", as.matrix(values))
		} else {
	    	values <- new("MxSymmetricSparse", matrix(values, nrow, ncol))
	    }
		retval <- callNextMethod(.Object, specification, values, name)
		return(retval)
	}
)

setMethod("omxVerifyMatrix", "SymmMatrix",
	function(.Object) {
		callNextMethod(.Object)
		verifySquare(.Object)
		values <- as.matrix(.Object@values)
		if (!all(values == t(values))) {
			stop(paste("Symmetric matrix",omxQuotes(.Object@name),"is not symmetric!"), 
				call.=FALSE)
		}
	}
)
