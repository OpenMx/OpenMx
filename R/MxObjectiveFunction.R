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
		F = "numeric"),
	contains = "MxObjective")

setMethod("initialize", "MxRAMObjective",
	function(.Object, A, S, F) {
		.Object@A <- A
		.Object@S <- S
		.Object@F <- F
		return(.Object)
	}
)

setClass(Class = "MxAlgebraObjective",
	representation = representation(
		algebra = "numeric"),
	contains = "MxObjective")

setMethod("initialize", "MxAlgebraObjective",
	function(.Object, algebra) {
		.Object@algebra <- algebra
		return(.Object)
	}
)

objectiveTypes <- c("FIML", "RAM", "algebra", "R")

mxFIMLObjective <- function(model, covariance, means) {
	if(class(model)[[1]] != "MxModel") {
		stop("First argument is not an MxModel object")
	}	
	if (typeof(covariance) != "character") {
		stop("Second argument is not a string (the name of the expected covariance matrix)")
	}
	if (typeof(means) != "character") {
		stop("Third argument is not a string (the name of the expected means matrix)")
	}
	covarianceMatrix <- match(covariance, names(model@matrices))
	meansMatrix <- match(means, names(model@matrices))
	if (is.na(covarianceMatrix)) {
		stop(paste("Could not find a matrix with name", covariance, "in the model."))
	}
	if (is.na(meansMatrix)) {
		stop(paste("Could not find a matrix with name", means, "in the model."))
	}
	return(new("MxFIMLObjective", - covarianceMatrix, - meansMatrix))
}

mxRAMObjective <- function(model, aMatrix = "A", sMatrix = "S", fMatrix = "F") {
	if(class(model)[[1]] != "MxModel") {
		stop("First argument is not an MxModel object")
	}	
	if (typeof(aMatrix) != "character") {
		stop("Second argument is not a string (the name of the 'A' matrix)")
	}	
	if (typeof(sMatrix) != "character") {
		stop("Third argument is not a string (the name of the 'S' matrix)")
	}
	if (typeof(fMatrix) != "character") {
		stop("Fourth argument is not a string (the name of the 'F' matrix)")
	}
	A <- match(aMatrix, names(model@matrices))
	S <- match(sMatrix, names(model@matrices))
	F <- match(fMatrix, names(model@matrices))	
	if (is.na(A)) {
		stop(paste("Could not find a matrix with name", aMatrix, "in the model."))
	}
	if (is.na(S)) {
		stop(paste("Could not find a matrix with name", sMatrix, "in the model."))
	}
	if (is.na(F)) {
		stop(paste("Could not find a matrix with name", fMatrix, "in the model."))
	}	
	return(new("MxRAMObjective", - A, - S, - F))
}

mxAlgebraObjective <- function(model, algebra) {
	if(class(model)[[1]] != "MxModel") {
		stop("First argument is not an MxModel object")
	}
	if (typeof(algebra) != "character") {
		stop("Second argument is not a string (the name of the algebra)")
	}
	algMatrix <- match(algebra, names(model@algebras))
	if (is.na(algMatrix)) {
		stop(paste("Could not find an algebra with name", algebra, "in the model."))
	}
	return(new("MxAlgebraObjective", algMatrix - 1))
}