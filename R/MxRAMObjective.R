setClass(Class = "MxRAMObjective",
	representation = representation(
		A = "MxCharOrNumber",
		S = "MxCharOrNumber",
		F = "MxCharOrNumber",
		data = "numeric"),
	contains = "MxBaseObjective")

setMethod("initialize", "MxRAMObjective",
	function(.Object, name, A, S, F, data = NA_real_) {
		.Object@name <- name
		.Object@A <- A
		.Object@S <- S
		.Object@F <- F
		.Object@data <- data
		return(.Object)
	}
)

setMethod("omxObjFunConvert", signature("MxRAMObjective", "MxFlatModel"), 
	function(.Object, model, definitions) {
		name <- .Object@name
		aMatrix <- .Object@A
		sMatrix <- .Object@S
		fMatrix <- .Object@F
		A <- omxLocateIndex(model, aMatrix, name)
		S <- omxLocateIndex(model, sMatrix, name)
		F <- omxLocateIndex(model, fMatrix, name)
		dIndex <- omxDataIndex(.Object@name, model@datasets)
		if (is.na(dIndex)) {
			msg <- paste("Could not find a data set for objective", 
				.Object@name, "in the model.")
			stop(msg)
		}
		return(new("MxRAMObjective", name, A, S, F, dIndex))
})

mxRAMObjective <- function(aMatrix = "A", sMatrix = "S", fMatrix = "F",
	name = omxUntitledName()) {
	if (typeof(name) != "character") {
		msg <- paste("Name argument is not a string",
			"(the name of the objective function)")
		stop(msg)
	}
	if (typeof(aMatrix) != "character") {
		msg <- paste("aMatrix argument is not a string",
			"(the name of the 'A' matrix)")
		stop(msg)
	}	
	if (typeof(sMatrix) != "character") {
		msg <- paste("sMatrix argument is not a string",
			"(the name of the 'S' matrix)")
		stop(msg)
	}
	if (typeof(fMatrix) != "character") {
		msg <- paste("fMatrix argument is not a string",
			"(the name of the 'F' matrix)")
		stop(msg)
	}
	return(new("MxRAMObjective", name, aMatrix, sMatrix, fMatrix))
}
