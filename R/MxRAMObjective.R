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
	A <- omxLocateIndex(model, aMatrix)
	S <- omxLocateIndex(model, sMatrix)
	F <- omxLocateIndex(model, fMatrix)
	if (is.na(A)) {
		stop(paste("Could not find a matrix/algebra with name", aMatrix, "in the model."))
	}
	if (is.na(S)) {
		stop(paste("Could not find a matrix/algebra with name", sMatrix, "in the model."))
	}
	if (is.na(F)) {
		stop(paste("Could not find a matrix/algebra with name", fMatrix, "in the model."))
	}	
	return(new("MxRAMObjective", A, S, F))
}