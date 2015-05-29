#' mxAvailableOptimizers
#'
#' List the Optimizers available in this version, e.g. "SLSQP" "CSOLNP"
#' 
#' note for advanced users: Special-purpose optimizers like Newton-Raphson or EM are not included in this list.
#'
#' @return - list of valid Optimizer names
#' @export
#' @family Miscellaneous Utility Functions
#' @seealso - \code{\link{mxOption}}(model, "Default optimizer")
#' @examples
#' mxAvailableOptimizers()
mxAvailableOptimizers <- function() {
	mxComputeGradientDescent()$availableEngines
}
