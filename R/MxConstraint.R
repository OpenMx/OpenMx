setClass(Class = "MxConstraint",
	representation = representation(
		name = "character",
		alg1 = "MxCharOrNumber",
		alg2 = "MxCharOrNumber",
		relation = "MxCharOrNumber"
	))
	
setMethod("initialize", "MxConstraint",
	function(.Object, name, alg1, alg2, relation) {
		.Object@name <- name
		.Object@alg1 <- alg1
		.Object@alg2 <- alg2
		.Object@relation <- relation
		return(.Object)
	}
)

mxConstraint <- function(name = omxUntitledName(),
	alg1, alg2, relation) {
	if (typeof(name) != "character") {
		stop(paste("Name argument is not a string",
		"(the name of the objective function)"))
	}
	if (missing(alg1) || typeof(alg1) != "character") {
		stop(paste("Alg1 argument is not a string",
		"(the name of the first algebra)"))		
	}		
	if (missing(alg2) || typeof(alg2) != "character") {
		stop(paste("Alg2 argument is not a string",
		"(the name of the second algebra)"))
	}
	if (missing(relation) || typeof(relation) != "character") {
		stop(paste("Relation argument is not a string",
		"(<, =, or >)"))
	}
	if (!(relation == "<" || 
		relation == "=" || relation == ">")) {
			stop("Relation must be <, =, or >")
		}
	return(new("MxConstraint", name, alg1, alg2, relation))
}