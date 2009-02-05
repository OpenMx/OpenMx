setClass(Class = "MxRAMObjective",
	representation = representation(
		A = "MxCharOrNumber",
		S = "MxCharOrNumber",
		F = "MxCharOrNumber"),
	contains = "MxBaseObjective")

setMethod("initialize", "MxRAMObjective",
	function(.Object, name, A, S, F) {
		.Object@name <- name
		.Object@A <- A
		.Object@S <- S
		.Object@F <- F
		return(.Object)
	}
)

setMethod("omxObjFunConvert", signature("MxRAMObjective", "MxModel"), 
	function(.Object, model) {
		name <- .Object@name
		aMatrix <- .Object@A
		sMatrix <- .Object@S
		fMatrix <- .Object@F
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
		return(new("MxRAMObjective", name, A, S, F))
})

mxRAMObjective <- function(name = omxUntitledName(), 
	aMatrix = "A", sMatrix = "S", fMatrix = "F") {
	if (typeof(name) != "character") {
		stop("Name argument is not a string (the name of the objective function)")
	}
	if (typeof(aMatrix) != "character") {
		stop("aMatrix argument is not a string (the name of the 'A' matrix)")
	}	
	if (typeof(sMatrix) != "character") {
		stop("sMatrix argument is not a string (the name of the 'S' matrix)")
	}
	if (typeof(fMatrix) != "character") {
		stop("fMatrix argument is not a string (the name of the 'F' matrix)")
	}
	return(new("MxRAMObjective", name, aMatrix, sMatrix, fMatrix))
}