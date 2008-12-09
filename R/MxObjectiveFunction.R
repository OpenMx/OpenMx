#
# Objective Functions are for optimizing!
#
setClass(Class = "MxObjective", 
	representation = representation("VIRTUAL"))

setClass(Class = "MxFIMLObjective",
	representation = representation(
		covariance = "SymmMatrix",
		means = "numeric"),
	contains = "MxObjective")

setMethod("initialize", "MxFIMLObjective",
	function(.Object, covariance, means) {
		.Object@covariance <- covariance
		.Object@means <- means
		return(.Object)
	}
)

setClass(Class = "MxCovObjective",
	representation = representation(
		covariance = "SymmMatrix"),
	contains = "MxObjective")

setMethod("initialize", "MxCovObjective",
	function(.Object, covariance) {
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

objectiveTypes <- c("FIML", "covariance", "algebra")

mxObjective <- function(type, covariance = NA, means = NA, algebra = NA) {
	if (is.na(match(type, objectiveTypes))) {
		stop(paste("Type must be one of:", paste(objectiveTypes, collapse=" ")))
	}
	if (type == "covariance") {
		if (single.na(covariance)) { stop("covariance must be specified") }
		return(new("MxCovObjective", covariance))
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
		

