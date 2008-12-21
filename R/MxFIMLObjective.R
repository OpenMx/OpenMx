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
	covarianceIndex <- omxLocateIndex(model, covariance)
	meansIndex <- omxLocateIndex(model, means)
	if (is.na(covarianceIndex)) {
		stop(paste("Could not find a matrix/algebra with name", covariance, "in the model."))
	}
	if (is.na(meansIndex)) {
		stop(paste("Could not find a matrix/algebra with name", means, "in the model."))
	}
	return(new("MxFIMLObjective", covarianceIndex, meansIndex))
}