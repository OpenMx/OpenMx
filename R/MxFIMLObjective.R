setClass(Class = "MxFIMLObjective",
	representation = representation(
		covariance = "MxCharOrNumber",
		means = "MxCharOrNumber",
		data = "numeric"),
	contains = "MxBaseObjective")

setMethod("initialize", "MxFIMLObjective",
	function(.Object, name, covariance, means, data = NA_real_) {
		.Object@name <- name
		.Object@covariance <- covariance
		.Object@means <- means
		.Object@data <- data
		return(.Object)
	}
)

setMethod("omxObjFunConvert", signature("MxFIMLObjective", "MxFlatModel"), 
	function(.Object, model) {
		name <- .Object@name
		covariance <- .Object@covariance
		means <- .Object@means
		covarianceIndex <- omxLocateIndex(model, covariance, name)
		meansIndex <- omxLocateIndex(model, means, name)
		dIndex <- omxDataIndex(.Object@name, model@datasets)
		if (is.na(covarianceIndex)) {
			stop(paste("Could not find a matrix/algebra with name", 
			covariance, "in the model."))
		}
		if (is.na(meansIndex)) {
			stop(paste("Could not find a matrix/algebra with name", 
			means, "in the model."))
		}
		if (is.na(dIndex)) {
			stop(paste("Could not find a data set for objective", 
			.Object@name, "in the model."))
		}
		return(new("MxFIMLObjective", name, covarianceIndex, meansIndex, dIndex))
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