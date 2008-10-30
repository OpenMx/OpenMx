setClassUnion("MxAlgebraFormula", c("call", "name", "logical"))

setClass(Class = "MxAlgebra",
	representation = representation(
		formula = "MxAlgebraFormula"))
		
setMethod("initialize", "MxAlgebra",
	function(.Object, formula = NA) {
		.Object@formula <- sys.call(which=-3)[[3]]
		return(.Object)
	}
)

mxAlgebra <- function(expression) {
	retval <- new("MxAlgebra", NA)
	retval@formula <- match.call()$expression
	return(retval)	
}

