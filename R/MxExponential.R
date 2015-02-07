##' Matrix exponential
##'
##' @param x matrix
##' @aliases omxExponential
expm <- function(x) {
	if (length(dim(x)) != 2 || nrow(x) != ncol(x)) {
		stop("input to omxExponential() must be a square matrix")
	}
	.Call(do_expm_eigen, x)
}

omxExponential <- expm

##' Matrix logarithm
##' 
##' @param x matrix
##' @param tol tolerance
logm <- function(x, tol = .Machine$double.eps) {
    d <- dim(x)
    if(length(d) != 2 || d[1] != d[2]) stop("'x' must be a quadratic matrix")
    .Call(do_logm_eigen, x)
}
