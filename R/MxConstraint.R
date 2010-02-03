#
#   Copyright 2007-2009 The OpenMx Project
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

mxConstraint <- function(alg1, relation, alg2, name = NA) {
	if (is.na(name)) {
		name <- omxUntitledName()
	}
	omxVerifyName(name)
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

omxConstraintRelations <- c("<", "=", ">")



displayMxConstraint <- function(object) {
	cat("MxConstraint", omxQuotes(object@name), '\n')
	cat("alg1 :", omxQuotes(object@alg1), '\n')
	cat("relation :", omxQuotes(object@relation), '\n')
	cat("alg2 :", omxQuotes(object@alg2), '\n')	
	invisible(object)
}

setMethod("print", "MxConstraint", function(x,...) { 
	displayMxConstraint(x) 
})

setMethod("show", "MxConstraint", function(object) { 
	displayMxConstraint(object) 
})