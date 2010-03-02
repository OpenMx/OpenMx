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

omxMnor <- function(cov, means, lbounds, ubounds) {
    
    cov <- as.matrix(cov)
    means <- as.matrix(means)
    lbounds <- as.matrix(lbounds)
    ubounds <- as.matrix(ubounds)

    if(nrow(cov) != ncol(cov)) {
        stop("Cov must be square")
    }
    if(nrow(means) > 1 && ncol(means) > 1) {
    	stop("'means' argument must be row or column vector")
    }
    if(nrow(lbounds) > 1 && ncol(lbounds) > 1) {
    	stop("'lbounds' argument must be row or column vector")    
    }
    if(nrow(ubounds) > 1 && ncol(ubounds) > 1) {
    	stop("'ubounds' argument must be row or column vector")    
    }
    
    if(ncol(cov) != length(means)) {
        stop("'means' must have length equal to diag(cov)")
    }
    if(ncol(cov) != length(lbounds)) {
        stop("'lbounds' must have length equal to diag(cov)")
    }
    if(ncol(cov) != length(ubounds)) {
        stop("'ubounds' must have length equal to diag(cov)")
    }
    
    retVal <- .Call("omxCallAlgebra", 
    	list(cov, means, lbounds, ubounds), 
    	omxLookupSymbolTable("omxMnor"), 
    	NA)
    return(as.matrix(as.numeric(retVal)))
    
}

omxAllInt <- function(cov, means, ...) {
    cov <- as.matrix(cov)
    means <- as.matrix(means)
    thresholdMats <- list(...)

    if(nrow(cov) != ncol(cov)) {
        stop("'cov' must be square")
    }
    if(nrow(means) > 1 && ncol(means) > 1) {
    	stop("'means' argument must be row or column vector")
    }
    if(ncol(cov) != length(means)) {
        stop("'means' must have length equal to diag(cov)")
    }
    if(sum(sapply(thresholdMats, ncol)) < ncol(cov)) {
        stop("'thresholds' must have at least as many total columns as 'cov'")
    }
    
    retVal <- .Call("omxCallAlgebra", 
        c(list(cov, means), thresholdMats),         # Flatten args into a single list
        omxLookupSymbolTable("omxAllInt"), 
        NA)
    
    return(as.matrix(as.numeric(retVal)))

}