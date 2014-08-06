setClass(Class = "MxFitFunctionMultigroup",
	 representation = representation(
	   groups = "MxOptionalCharOrNumber",
	     verbose= "integer"),
	 contains = "MxBaseFitFunction")

setMethod("initialize", "MxFitFunctionMultigroup",
	function(.Object, groups, verbose, name = 'fitfunction') {
		.Object@name <- name
		.Object@groups <- groups
		.Object@verbose <- verbose
		return(.Object)
	}
)

setMethod("genericFitDependencies", signature("MxFitFunctionMultigroup"),
	function(.Object, flatModel, dependencies) {
	dependencies <- callNextMethod()
	dependencies <- imxAddDependency(.Object@groups, .Object@name, dependencies)
	return(dependencies)
})

setMethod("qualifyNames", signature("MxFitFunctionMultigroup"), 
	function(.Object, modelname, namespace) {
		.Object@name <- imxIdentifier(modelname, .Object@name)
		return(.Object)
})

# "model.algebra" or "model" for "model.fitfunction"
setMethod("genericFitFunConvert", "MxFitFunctionMultigroup", 
	function(.Object, flatModel, model, labelsData, defVars, dependencies) {
		name <- .Object@name
		if (length(.Object@groups)) .Object@groups <- vapply(.Object@groups, function(group) {
			path <- unlist(strsplit(group, imxSeparatorChar, fixed = TRUE))
			if (length(path) == 1) {
				ff <- paste(path, "fitfunction", sep=".")
				length(model@algebras) + imxLocateIndex(flatModel, ff, name)
			} else if (length(path) == 2) {
				# restrict to algebra or fitfunction TODO
				imxLocateIndex(flatModel, group, name)
			}
		}, 1L)
		return(.Object)
})

setMethod("generateReferenceModels", "MxFitFunctionMultigroup",
	function(.Object, model) {
		grpnames <- unlist(strsplit(model$fitfunction$groups, split=".fitfunction", fixed=TRUE))
		grpmodels <- list()
		for(i in 1:length(grpnames)){
			grpmodels[[i]] <- ReferenceModelHelper(model[[ grpnames[i] ]], run=FALSE)
		}
		sgrpmodels <- sapply(grpmodels, "[[", 1) #extract saturated models
		sgrpfits <- mxFitFunctionMultigroup(paste(sapply(sgrpmodels, slot, name="name"), ".fitfunction", sep=""))
		saturatedModel <- mxModel(name=paste("Saturated", modelName), sgrpmodels, sgrpfits)
		igrpmodels <- sapply(grpmodels, "[[", 2) #extract independence models
		igrpfits <- mxFitFunctionMultigroup(paste(sapply(igrpmodels, slot, name="name"), ".fitfunction", sep=""))
		independenceModel <- mxModel(name=paste("Independence", modelName), igrpmodels, igrpfits)
		return(list(Saturated=saturatedModel, Independence=independenceModel))
	})

mxFitFunctionMultigroup <- function(groups, ..., verbose=0L) {
	garbageArguments <- list(...)
	if (length(garbageArguments) > 0) {
		stop("mxFitFunctionMultigroup does not accept values for the '...' argument")
	}

	if (length(groups) == 0) stop("mxFitFunctionMultigroup: at least 1 fitfunction must be provided")

	return(new("MxFitFunctionMultigroup", groups, verbose))
}
