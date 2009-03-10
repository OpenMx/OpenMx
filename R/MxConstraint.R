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

mxConstraint <- function(alg1, alg2, relation,
	name = omxUntitledName()) {
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
	if (!(relation %in% omxConstraintRelations)) {
			clist <- paste(omxConstraintRelations, 
				collapse = ", ")
			msg <- paste("Relation must be in the list:",
				clist)
			stop(msg)
		}
	return(new("MxConstraint", name, alg1, alg2, relation))
}

omxConvertConstraints <- function(flatModel) {
	return(lapply(flatModel@constraints, function(x) {
		omxConvertSingleConstraint(x, flatModel)}))
}


omxConvertSingleConstraint <- function(constraint, flatModel) {
	index1 <- omxLocateIndex(flatModel, constraint@alg1,
		constraint@name)
	index2 <- omxLocateIndex(flatModel, constraint@alg2,
		constraint@name)
	index3 <- match(constraint@relation, omxConstraintRelations)
	if(is.na(index3)) {
		clist <- paste(omxConstraintRelations, 
				collapse = ",")		
		msg <- paste("The relation for constraint", 
			omxQuotes(constraint@name),
			"is not in the following list:",
			clist)
		stop(msg, call.=FALSE)
	}
	return(list(index1,index2,index3 - 1))
}

omxConstraintRelations <- c("<", "=", ">")   
