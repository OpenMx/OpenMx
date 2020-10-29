setClass(Class = "MxExpectationMixture",
	representation = representation(
		weights = "MxCharOrNumber",
		components = "MxOptionalCharOrNumber",
		verbose = "integer",
		scale = "character"
	),
	contains = "MxBaseExpectation")

setMethod("initialize", "MxExpectationMixture",
          function(.Object, components, weights, verbose, scale, name = 'expectation') {
		  .Object@data = as.integer(NA)
		  .Object@name <- name
		  .Object@components <- components
		  .Object@weights <- weights
		  .Object@verbose <- verbose
		  .Object@scale <- scale
		  .Object
	  })

setMethod("genericExpDependencies", signature("MxExpectationMixture"),
	  function(.Object, dependencies) {
		  components <- paste(.Object@components, "expectation", sep=".")
		  sources <- c(.Object@weights, components)
		  dependencies <- imxAddDependency(sources, .Object@name, dependencies)
		  return(dependencies)
	  })

setMethod("genericExpFunConvert", signature("MxExpectationMixture"), 
	  function(.Object, flatModel, model, labelsData, dependencies) {
		  if (length(.Object@components)) {
			  origComponents <- .Object@components
			  .Object@components <- vapply(.Object@components, function(group) {
				  eNumber <- match(paste(group, "expectation", sep="."),
						   names(flatModel@expectations))
				  eNumber - 1L
			  }, 1L, USE.NAMES = FALSE)
			  if (any(is.na(.Object@components))) {
				  stop(paste(model@name,": cannot locate expectation ",
					     omxQuotes(origComponents[is.na(.Object@components)]), sep=""),
				       call. = FALSE)
			  }
		  }
		  .Object
	  })

setMethod("qualifyNames", signature("MxExpectationMixture"), 
	function(.Object, modelname, namespace) {
		.Object@name <- imxIdentifier(modelname, .Object@name)
		.Object@data <- imxConvertIdentifier(.Object@data, 
			modelname, namespace)
		for (s in c('weights')) {
			if (is.null(slot(.Object, s))) next;
			slot(.Object, s) <-
				imxConvertIdentifier(slot(.Object, s), modelname, namespace)
		}
		.Object
})

setMethod("genericNameToNumber", signature("MxExpectationMixture"),
	  function(.Object, flatModel, model) {
		  name <- .Object@name
		  .Object@data <- imxLocateIndex(flatModel, .Object@data, name)
		  .Object@weights <- imxLocateIndex(flatModel, .Object@weights, name)
		  .Object
	  })

setMethod("genericGetExpected", signature("MxExpectationMixture"),
	function(.Object, model, what, defvar.row=1, subname=model@name) {
		ret <- list()
		if ('weights' %in% what) {
			wname <- .modifyDottedName(model@name, .Object@weights)
			weights <- mxEvalByName(wname, model, compute=TRUE, defvar.row=defvar.row)
			if (.Object@scale == 'softmax') {
				weights <- exp(weights)
			} else if (.Object@scale == 'sum') {
				#OK
			} else if(.Object@scale == 'none') {
				#OK
			} else { stop(.Object@scale) }
			if(.Object@scale != 'none'){
				weights <- weights / rowSums(weights)
			}
			ret[['weights']] <- weights
		}
		ret
	})

setMethod("genericGenerateData", signature("MxExpectationMixture"),
	function(.Object, model, nrows, subname, empirical, returnModel, use.miss,
		   .backend, nrowsProportion) {
		origData <- findDataForSubmodel(model, subname)
		origRows <- if (!is.null(origData)) { nrowMxData(origData) } else { NULL }
		nrows <- calcNumRows(nrows, nrowsProportion, origRows, subname)

		cdata <- list()
		for (c1 in .Object@components) {
		  cdata[[c1]] <- mxGenerateData(model, returnModel=FALSE, nrows=nrows, use.miss=FALSE,
						subname=c1, empirical=empirical, .backend=.backend)
		}
		data <- cdata[[1]]

		# This is an inefficient way to generate data. It would be
		# better to generate cpick first and then generate 1 row
		# at a time from the component expectations. I didn't code
		# it that way because the API is not really set up for
		# generating data 1 row at a time.
		cpick <- NULL
		doDefVar <- imxHasDefinitionVariable(model)
		if(doDefVar){
			if (origData$type != 'raw') {
				stop(paste("Definition variable(s) found, but original data is type",
					omxQuotes(origData$type)))
			}
			origData <- origData$observed
			if(nrows != nrow(origData)){
				stop("Definition variable(s) found, but the number of rows in the data do not match the number of rows requested for data generation.")
			}
			cpick <- rep(NA, nrows)
			for (rx in 1:nrows) {
				weights <- mxGetExpected(model, "weights", defvar.row=rx)
				cpick[rx] <- sample.int(length(.Object@components), 1, prob=weights)
			}
		} else {
			weights <- mxGetExpected(model, "weights")
			cpick <- sample.int(length(.Object@components), nrows, replace=TRUE, prob=weights)
		}
		if (length(.Object@components) > 1) for (cx in 2:length(.Object@components)) {
			data[cpick==cx,] <- cdata[[cx]][cpick == cx,]
		}
		if(doDefVar){
			for (dcol in setdiff(colnames(origData), colnames(data))) {
				data[[dcol]] <- origData[[dcol]]
			}
		}
		if (returnModel) {
		  mxModel(model[[subname]], mxData(as.data.frame(data), "raw"))
		} else {
		  as.data.frame(data)
		}
	})

mxExpectationMixture <- function(components, weights="weights",
				      ..., verbose=0L, scale=c('softmax', 'sum', 'none')) {
  prohibitDotdotdot(list(...))

	scale <- match.arg(scale)

	new("MxExpectationMixture", components, weights,
	    as.integer(verbose), scale)
}
