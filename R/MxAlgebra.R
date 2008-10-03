setClass(Class = "MxAlgebra",
	representation = representation(
		formula = "expression"))
		
setMethod("initialize", "MxAlgebra",
	function(.Object, formula = NA) {
		.Object@formula <- as.expression(sys.call(which=-3)[[3]])
		return(.Object)
	}
)