setClass(Class = "MxFitFunctionMultigroup",
	 representation = representation(
	   groups = "MxOptionalCharOrNumber",
	     verbose= "integer"),
	 contains = "MxBaseFitFunction")

setMethod("initialize", "MxFitFunctionMultigroup",
	function(.Object, ...) {
    .Object <- callNextMethod()
		.Object@groups <- ..1
		.Object@verbose <- ..2
		.Object
	})

setMethod("genericFitDependencies", signature("MxFitFunctionMultigroup"),
	function(.Object, flatModel, dependencies) {
	dependencies <- callNextMethod()
	groups <- vapply(.Object@groups, function(group) {
		path <- unlist(strsplit(group, imxSeparatorChar, fixed = TRUE))
		if (length(path) == 1) {
			group <- paste(path, "fitfunction", sep=".")
		}
		group
	}, "")
	dependencies <- imxAddDependency(groups, .Object@name, dependencies)
	return(dependencies)
})

setMethod("genericFitRename", signature("MxFitFunctionMultigroup"),
	function(.Object, oldname, newname) {
		.Object@groups <- sapply(.Object@groups, renameUnqualifiedReference, oldname, newname)
		.Object
})

# "model.algebra" or "model" for "model.fitfunction"
setMethod("genericFitFunConvert", "MxFitFunctionMultigroup",
	function(.Object, flatModel, model, labelsData, dependencies) {
    .Object <- callNextMethod()
		name <- .Object@name
		if (length(.Object@groups)) {
			origGroups <- .Object@groups
			.Object@groups <- vapply(.Object@groups, function(group) {
				path <- unlist(strsplit(group, imxSeparatorChar, fixed = TRUE))
				if (length(path) == 1) {
					group <- paste(path, "fitfunction", sep=".")
				}
				algebraNumber <- match(group, append(names(flatModel@algebras),
								     names(flatModel@fitfunctions)))
				algebraNumber - 1L
			}, 1L)
			if (any(is.na(.Object@groups))) {
				stop(paste(name,": cannot locate algebra/fitfunction ",
					   omxQuotes(origGroups[is.na(.Object@groups)]), sep=""),
				     call. = FALSE)
			}
		}
		return(.Object)
})

aggregateSubrefmodels <- function(modelName, submodels) {
	if (is(submodels[[1]], "MxModel")) {
		fit <- mxFitFunctionMultigroup(paste(sapply(submodels, slot, name="name"), ".fitfunction", sep=""))
		mxModel(name=modelName, submodels, fit)
	} else if (is.numeric(submodels[[1]])) {
		list(sum(unlist(submodels[1,])),
		     sum(unlist(submodels[2,])))
	} else {
		stop(paste("Not sure how to aggregate:\n", paste(deparse(submodels), collapse="\n")))
	}
}

setMethod("generateReferenceModels", "MxFitFunctionMultigroup",
	function(.Object, model, distribution, equateThresholds) {
		grpnames <- unlist(strsplit(model$fitfunction$groups, split=".fitfunction", fixed=TRUE))
		grpmodels <- list()
		for(i in 1:length(grpnames)){
			grpmodels[[i]] <- ReferenceModelHelper(model[[ grpnames[i] ]], distribution, equateThresholds)
		}
		sgrpmodels <- sapply(grpmodels, "[[", 1)
		saturatedModel <- aggregateSubrefmodels(paste("Saturated", model@name), sgrpmodels)
		igrpmodels <- sapply(grpmodels, "[[", 2)
		independenceModel <- aggregateSubrefmodels(paste("Independence", model@name), igrpmodels)
		return(list(Saturated=saturatedModel, Independence=independenceModel))
	})

mxFitFunctionMultigroup <- function(groups, ..., verbose=0L) {
	prohibitDotdotdot(list(...))

	if (length(groups) == 0) stop("mxFitFunctionMultigroup: at least 1 fitfunction must be provided")

	return(new("MxFitFunctionMultigroup", groups, as.integer(verbose)))
}
