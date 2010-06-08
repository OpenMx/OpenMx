#
#   Copyright 2007-2010 The OpenMx Project
#
#   Licensed under the Apache License, Version 2.0 (the "License");
#   you may not use this file except in compliance with the License.
#   You may obtain a copy of the License at
# 
#        http://www.apache.org/licenses/LICENSE-2.0
# 
#    Unless required by applicable law or agreed to in writing, software
#    distributed under the License is distributed on an "AS IS" BASIS,
#    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#   See the License for the specific language governing permissions and
#   limitations under the License.


setClass(Class = "MxConstraint",
	representation = representation(
		name = "character",
		formula = "MxAlgebraFormula",
		alg1 = "MxCharOrNumber",
		alg2 = "MxCharOrNumber",
		relation = "MxCharOrNumber"
	))
	
setMethod("initialize", "MxConstraint",
	function(.Object, name, formula) {
		.Object@name <- name
		.Object@formula <- formula
		.Object@alg1 <- as.character(NA)
		.Object@alg2 <- as.character(NA)
		.Object@relation <- as.character(NA)
		return(.Object)
	}
)

mxConstraint <- function(expression, name = NA) {
	if (is.na(name)) {
		name <- omxUntitledName()
	}
	omxVerifyName(name, 0)
	if (typeof(name) != "character") {
		stop(paste("Name argument is not a string",
		"(the name of the objective function)"))
	}
	formula <- match.call()$expression
	if (length(formula) != 3) {
		stop(paste("A constraint must be of the form:",
			omxQuotes("exp1 [<, ==, >] exp2")),
			call. = FALSE)
	}
	relation <- as.character(formula[[1]])
	if (!(relation %in% omxConstraintRelations)) {
		stop(paste("A constraint must be of the form:",
			omxQuotes("exp1 [<, ==, >] exp2")),
			call. = FALSE)
	}
	return(new("MxConstraint", name, formula))
}

convertConstraints <- function(flatModel) {
	return(lapply(flatModel@constraints, function(x) {
		convertSingleConstraint(x, flatModel)}))
}


convertSingleConstraint <- function(constraint, flatModel) {
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

omxConstraintRelations <- c("<", "==", ">")



displayConstraint <- function(object) {
	cat("MxConstraint", omxQuotes(object@name), '\n')
	cat("@formula: ", deparse(object@formula, width.cutoff=500L), '\n')
}

setMethod("print", "MxConstraint", function(x,...) { displayConstraint(x) })
setMethod("show", "MxConstraint", function(object) { displayConstraint(object) })
