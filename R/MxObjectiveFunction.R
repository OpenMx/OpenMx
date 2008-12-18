#
# Objective Functions are for optimizing!
#
setClass(Class = "MxObjective", 
	representation = representation("VIRTUAL"))

setClass(Class = "MxFIMLObjective",
	representation = representation(
		covariance = "numeric",
		means = "numeric"),
	contains = "MxObjective")

setMethod("initialize", "MxFIMLObjective",
	function(.Object, covariance, means) {
		.Object@covariance <- covariance
		.Object@means <- means
		return(.Object)
	}
)

setClass(Class = "MxRAMObjective",
	representation = representation(
		A = "numeric",
		S = "numeric",
		F = "numeric",
		covariance = "SymmMatrix"),
	contains = "MxObjective")

setMethod("initialize", "MxRAMObjective",
	function(.Object, A = -1, S = -2, F = -3, covariance) {
		.Object@A <- A
		.Object@S <- S
		.Object@F <- F
		.Object@covariance <- covariance
		return(.Object)
	}
)

setClass(Class = "MxAlgebraObjective",
	representation = representation(
		algebra = "MxAlgebra"),
	contains = "MxObjective")

setMethod("initialize", "MxAlgebraObjective",
	function(.Object, algebra) {
		.Object@algebra <- algebra
		return(.Object)
	}
)

objectiveTypes <- c("FIML", "RAM", "algebra", "R")

mxObjective <- function(type, covariance = NA, means = NA, algebra = NA) {
	if (is.na(match(type, objectiveTypes))) {
		stop(paste("Type must be one of:", paste(objectiveTypes, collapse=" ")))
	}
	if (type == "RAM") {
		if (single.na(covariance)) { stop("covariance must be specified") }
		return(new("MxRAMObjective", -1, -2, -3, covariance))
	} else if (type == "FIML") {
		if (single.na(means)) { stop("means must be specified") }
		if (single.na(covariance)) { stop("covariance must be specified") }
		return(new("MxFIMLObjective", covariance, means))
	} else if (type == "algebra") {
		if (!is(algebra, "MxAlgebra")) {
			stop("algebra must be of type MxAlgebra")
		}
		return(new("MxAlgebraObjective", algebra))
	}
}


