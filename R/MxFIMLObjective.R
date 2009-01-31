setClass(Class = "MxFIMLObjective",
	representation = representation(
		covariance = "MxCharOrNumber",
		means = "MxCharOrNumber"),
	contains = "MxObjective")

setMethod("initialize", "MxFIMLObjective",
	function(.Object, name, covariance, means) {
		.Object@name <- name
		.Object@covariance <- covariance
		.Object@means <- means
		return(.Object)
	}
)

setMethod("omxObjFunConvert", signature("MxFIMLObjective", "MxModel"), 
	function(.Object, model) {
		name <- .Object@name
		covariance <- .Object@covariance
		means <- .Object@means
		covarianceIndex <- omxLocateIndex(model, covariance)
		meansIndex <- omxLocateIndex(model, means)
		if (is.na(covarianceIndex)) {
			stop(paste("Could not find a matrix/algebra with name", 
			covariance, "in the model."))
		}
		if (is.na(meansIndex)) {
			stop(paste("Could not find a matrix/algebra with name", 
			means, "in the model."))
		}
		return(new("MxFIMLObjective", name, covarianceIndex, meansIndex))
})

mxFIMLObjective <- function(name = omxUntitledName(), covariance, means) {
	if (typeof(name) != "character") {
		stop("Name argument is not a string (the name of the objective function)")
	}
	if (missing(covariance) || typeof(covariance) != "character") {
		stop("Covariance argument is not a string (the name of the expected covariance matrix)")
	}
	if (missing(means) || typeof(means) != "character") {
		stop("Means argument is not a string (the name of the expected means matrix)")
	}
	return(new("MxFIMLObjective", name, covariance, means))
}