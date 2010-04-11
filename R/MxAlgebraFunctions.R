#
#   Copyright 2007-2010 The OpenMx Project
#
#   Licensed under the Apache License, Version 2.0 (the "License");
#   you may not use this file except in compliance with the License.
#   You may obtain a copy of the License at
# 
#        http://www.apache.org/licenses/LICENSE-2.0
# 
#   Unless required by applicable law or agreed to in writing, software
#   distributed under the License is distributed on an "AS IS" BASIS,
#   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#   See the License for the specific language governing permissions and
#   limitations under the License.


tr <- function(x) {
	if (is.matrix(x)) {
		return(sum(diag(x)))
	} else {
		return(NA)
	}
}

"%&%" <- function(x, y) {
	return(x %*% y %*% t(x))
}

"%^%" <- function(x, y) {
	return(kronecker(x, y, '^'))
}

cvectorize <- function(x) {	
	return(matrix(x, length(x), 1))
}

rvectorize <- function(x) {
	return(matrix(t(x), length(x), 1))
}

vech <- function(x) {
	return(x[lower.tri(x, diag=TRUE)])
}

vechs <- function(x) {
	return(x[lower.tri(x, diag=FALSE)])
}

diag2vec <- function(x) {
	return(as.matrix(diag(as.matrix(x))))
}

vec2diag <- function(x) {
	x <- as.matrix(x)
	if (nrow(x) != 1 && ncol(x) != 1) {
		stop("argument must be a row or column vector")
	}
	if (nrow(x) * ncol(x) == 1) {
		return(x)
	} else {
		return(as.matrix(diag(as.numeric(x))))
	}
}

omxLookupSymbolTable <- function(name) {
	index <- which(omxSymbolTable["R.name"] == name)
	if(length(index) == 0) {
		stop(paste("Internal error, function",
			name, "cannot be found in OpenMx symbol table"),
			call. = FALSE)
	} else if (length(index) > 1) {
		stop(paste("Internal error, function",
			name, "appears twice in OpenMx symbol table"),
			call. = FALSE)	
	}
	return(as.integer(index - 1))
}

omxMnor <- function(covariance, means, lbound, ubound) {
    
    covariance <- as.matrix(covariance)
    means <- as.matrix(means)
    lbound <- as.matrix(lbound)
    ubound <- as.matrix(ubound)

    if(nrow(covariance) != ncol(covariance)) {
        stop("Covariance must be square")
    }
    if(nrow(means) > 1 && ncol(means) > 1) {
    	stop("'means' argument must be row or column vector")
    }
    if(nrow(lbound) > 1 && ncol(lbound) > 1) {
    	stop("'lbound' argument must be row or column vector")    
    }
    if(nrow(ubound) > 1 && ncol(ubound) > 1) {
    	stop("'ubound' argument must be row or column vector")    
    }
    
    if(ncol(covariance) != length(means)) {
        stop("'means' must have length equal to diag(covariance)")
    }
    if(ncol(covariance) != length(lbound)) {
        stop("'lbound' must have length equal to diag(covariance)")
    }
    if(ncol(covariance) != length(ubound)) {
        stop("'ubound' must have length equal to diag(covariance)")
    }
    
    retVal <- .Call("omxCallAlgebra", 
    	list(covariance, means, lbound, ubound), 
    	omxLookupSymbolTable("omxMnor"), 
    	NA)
    return(as.matrix(as.numeric(retVal)))
    
}

omxAllInt <- function(covariance, means, ...) {
    covariance <- as.matrix(covariance)
    means <- as.matrix(means)
    thresholdMats <- list(...)

    if(nrow(covariance) != ncol(covariance)) {
        stop("'covariance' must be square")
    }
    if(nrow(means) > 1 && ncol(means) > 1) {
    	stop("'means' argument must be row or column vector")
    }
    if(ncol(covariance) != length(means)) {
        stop("'means' must have length equal to diag(cov)")
    }
    
    if(sum(sapply(thresholdMats, ncol)) < ncol(covariance)) {
        stop("'thresholds' must have at least as many total columns as 'covariance'")
    }

    if(min(sapply(thresholdMats, nrow)) < 2) {
        stop("every column of 'thresholds' must have at least two rows: one lower bound and one upper")
    }
    
    retVal <- .Call("omxCallAlgebra", 
        c(list(covariance, means), thresholdMats),         # Flatten args into a single list
        omxLookupSymbolTable("omxAllInt"), 
        NA)
    
    return(as.matrix(as.numeric(retVal)))

}

eigenvec <- function(x) {
    x <- as.matrix(x)
    if(nrow(x) != ncol(x)) {
        stop("matrix must be square")
    }
    
    retval <- .Call("omxCallAlgebra", 
        list(x),         # Flatten args into a single list
        omxLookupSymbolTable("eigenvec"), 
        NA)
        
    return(matrix(as.numeric(retval), nrow(x), ncol(x)))
}

ieigenvec <- function(x) {
    x <- as.matrix(x)
    if(nrow(x) != ncol(x)) {
        stop("matrix must be square")
    }
    
    retval <- .Call("omxCallAlgebra", 
        list(x),         # Flatten args into a single list
        omxLookupSymbolTable("ieigenvec"), 
        NA)
        
    return(matrix(as.numeric(retval), nrow(x), ncol(x)))
}

eigenval <- function(x) {
    x <- as.matrix(x)
    if(nrow(x) != ncol(x)) {
        stop("matrix must be square")
    }
    
    retval <- .Call("omxCallAlgebra", 
        list(x),         # Flatten args into a single list
        omxLookupSymbolTable("eigenval"), 
        NA)

    return(as.matrix(as.numeric(retval)))
}

ieigenval <- function(x) {
    x <- as.matrix(x)
    if(nrow(x) != ncol(x)) {
        stop("matrix must be square")
    }
    
    retval <- .Call("omxCallAlgebra", 
        list(x),         # Flatten args into a single list
        omxLookupSymbolTable("ieigenval"), 
        NA)
        
    return(as.matrix(as.numeric(retval)))
}
