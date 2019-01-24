#
#   Copyright 2007-2019 by the individuals mentioned in the source code history
#
#   Licensed under the Apache License, Version 2.0 (the "License");
#   you may not use this file except in compliance with the License.
#   You may obtain a copy of the License at
# 
#        http://www.apache.org/licenses/LICENSE-2.0
# 
#   Unless required by applicable law or agreed to in writing, software
#   distributed under the License is distributed on an "AS IS" BASIS,
#   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#   See the License for the specific language governing permissions and
#   limitations under the License.

##' @name MxDataFrameOrMatrix-class
##' @rdname MxDataFrameOrMatrix-class
##' @title MxDataFrameOrMatrix
##'
##' @description
##' Internal class that is the union of data.frame and matrix.
##'
##' @aliases
##' MxDataFrameOrMatrix
##' @details
##' Not to be used.
setClassUnion("MxDataFrameOrMatrix", c("data.frame", "matrix"))

setClass(Class = "NonNullData",
	representation = representation(
		verbose = "integer"))

##' @name MxDataStatic-class
##' @rdname MxDataStatic-class
##' @title Create static data
##'
##' @description
##' Internal static data class.
##'
##' @aliases
##' MxDataStatic
##' print,MxDataStatic-method
##' show,MxDataStatic-method
##' @details
##' Not to be used.
setClass(Class = "MxDataStatic",
	 contains = "NonNullData",
	 representation = representation(
		observed = "MxDataFrameOrMatrix",
		means  = "matrix",
		type   = "character",
		numObs = "numeric",
		observedStats = "list",
		.isSorted = "logical",  # remove slot TODO
		.needSort = "logical",
		.parallel = "logical",
	     primaryKey = "MxCharOrNumber",
	     weight = "MxCharOrNumber",
	     frequency = "MxCharOrNumber",
		name   = "character"))

setClass(Class = "MxDataDynamic",
	 contains = "NonNullData",
	 representation = representation(
	     type        = "character",
	     expectation = "MxCharOrNumber",
	     numObs = "numeric",             # output
	     name        = "character"))

setClassUnion("MxData", c("NULL", "MxDataStatic", "MxDataDynamic"))

setMethod("initialize", "MxDataStatic",
	  function(.Object, observed, means, type, numObs, observedStats,
		   sort, primaryKey, weight, frequency, verbose, .parallel) {
		.Object@observed <- observed
		.Object@means <- means
		.Object@type <- type
		.Object@numObs <- numObs
		.Object@observedStats <- observedStats
		.Object@name <- "data"
		if (is.na(primaryKey)) {
			if (is.na(sort)) {
				sort <- TRUE
			}
		} else {
			if (!is.na(sort) && sort) {
				warning("sort=TRUE is not supported when a primary key is provided")
			}
			sort <- FALSE
		}
		.Object@.isSorted <- FALSE
		.Object@.needSort <- sort
		.Object@.parallel <- .parallel
		.Object@primaryKey <- primaryKey
		.Object@weight <- weight
		.Object@frequency <- frequency
		.Object@verbose <- verbose
		return(.Object)
	}
)

setMethod("initialize", "MxDataDynamic",
	function(.Object, type, expectation, verbose, name = "data") {
		.Object@type <- type
		.Object@expectation <- expectation
		.Object@verbose <- verbose
		.Object@name <- name
		return(.Object)
	}
)

setMethod("$", "MxData", imxExtractSlot)

setReplaceMethod("$", "MxData",
       function(x, name, value) {
               return(imxReplaceSlot(x, name, value))
       }
)

setMethod("names", "MxData", slotNames)

##' Valid types of data that can be contained by MxData
imxDataTypes <- c("raw", "cov", "cor", "sscp")

##' Create dynamic data
##'
##' @param type type of data
##' @param ...  Not used.  Forces remaining arguments to be specified by name.
##' @param expectation the name of the expectation to provide the data
##' @param verbose Increase runtime debugging output
##' @aliases
##' MxDataDynamic-class
##' print,MxDataDynamic-method
##' show,MxDataDynamic-method
mxDataDynamic <- function(type, ..., expectation, verbose=0L) {
	garbageArguments <- list(...)
	if (length(garbageArguments) > 0) {
		stop("mxDataDynamic does not accept values for the '...' argument")
	}
	if (type != "cov") stop("Type must be set to 'cov'")
	verbose <- as.integer(verbose)
	return(new("MxDataDynamic", type, expectation, verbose))
}

mxData <- function(observed, type, means = NA, numObs = NA, acov=NA, fullWeight=NA,
		   thresholds=NA, ...,
		   observedStats=NA, sort=NA, primaryKey = as.character(NA), weight = as.character(NA),
		   frequency = as.character(NA), verbose=0L, .parallel=TRUE) {
	garbageArguments <- list(...)
	if (length(garbageArguments) > 0) {
		stop("mxData does not accept values for the '...' argument")
	}
	if (length(means) == 1 && is.na(means)) means <- as.numeric(NA)
	if (missing(observedStats)) {
		observedStats <- list()
		if (!missing(acov)) observedStats <- c(observedStats, acov=acov)
		if (!missing(fullWeight)) observedStats <- c(observedStats, fullWeight=fullWeight)
		if (!missing(thresholds)) observedStats <- c(observedStats, thresholds=thresholds)
	} else {
		if (!missing(acov) || !missing(fullWeight) || !missing(thresholds)) {
			stop("acov, fullWeight, and thresholds must be passed in the observedStats list")
		}
	}
	rm(acov)
	rm(fullWeight)
	rm(thresholds)
	if (missing(observed) || !is(observed, "MxDataFrameOrMatrix")) {
		stop("Observed argument is neither a data frame nor a matrix")
	}
	dups <- duplicated(colnames(observed))
	if (any(dups)) {
		stop(paste("Column names must be unique. Duplicated:",
			   omxQuotes(colnames(observed)[dups])))
	}
	if (missing(type) || (!is.character(type)) || (length(type) > 1) || 
		is.na(match(type, imxDataTypes))) {
			numTypes = length(imxDataTypes)
			stop(paste("Type must be set to one of: ", "'", 
			paste(imxDataTypes[1:(numTypes-1)], collapse="' '"),
			"' or '", imxDataTypes[numTypes], "'", sep=""))
	}
	if (type == "sscp") {
		stop(paste("'sscp' is not yet implemented."))
	}
	if ((!is.vector(means) && !(prod(dim(means)) == length(means))) || !is.numeric(means)) {
		stop("Means argument must be of numeric vector type")
	}
	if (type != "raw" && is.na(numObs)) {
		stop("Number of observations must be specified for non-raw data, i.e., add numObs=XXX to mxData()")
	}
	if (type == "cov") {
		verifyCovarianceMatrix(observed)
	}
	if (type == "cor") {
		verifyCorrelationMatrix(observed)
	}
	if (!is.null(observedStats[['cov']])) {
		verifyCovarianceMatrix(observedStats[['cov']], nameMatrix="covariance")
	}
	if (!is.null(observedStats[['acov']])) {
		verifyCovarianceMatrix(observedStats[['acov']], nameMatrix="asymptotic")
	}
	fw <- observedStats[['fullWeight']]
	if(!is.null(fw)) {
		verifyCovarianceMatrix(fw, nameMatrix="asymptotic")
	}
	thr <- observedStats[['thresholds']]
	if (!is.null(thr)) {
		verifyThresholdNames(thr, observed)
	}
	lapply(dimnames(observed)[[2]], imxVerifyName, -1)
	if(is.matrix(means)){meanNames <- colnames(means)} else {meanNames <- names(means)}
	means <- as.matrix(means)
	dim(means) <- c(1, length(means))
	colnames(means) <- meanNames

	if (length(primaryKey) > 1) {
		stop("More than 1 primary key is not implemented yet")
	} else if (!is.na(primaryKey)) {
		if (type != "raw") {
			stop(paste("Raw data is required when a primary key is provided"))
		}
		if (!(primaryKey %in% colnames(observed))) {
			stop(paste("Primary key", omxQuotes(primaryKey),
				   "must be provided in the observed data"))
		}
	}

	if (length(weight) > 1) {
		stop("Only one weight column can be specified")
	} else if (!is.na(weight)) {
		if (type != "raw") {
			stop(paste("Raw data is required when a weight column is provided"))
		}
		if (!(weight %in% colnames(observed))) {
			stop(paste("Weight column", omxQuotes(weight),
				   "must be provided in the observed data"))
		}
	}
	if (length(frequency) > 1) {
		stop("Only one frequency column can be specified")
	} else if (!is.na(frequency)) {
		if (type != "raw") {
			stop(paste("Raw data is required when a frequency column is provided"))
		}
		if (!(frequency %in% colnames(observed))) {
			stop(paste("Frequency column", omxQuotes(frequency),
				   "must be provided in the observed data"))
		}
	}
	if (type == "raw") {
		if (!is.na(frequency)) {
			obsCount <- sum(observed[,frequency])
		} else {
			obsCount <- nrow(observed)
		}

		if (missing(numObs)) {
			numObs <- obsCount
		} else {
			if (numObs != obsCount) {
				# stop?
				warning(paste("numObs of", numObs, "does not match the number of observations",
					      "found in the observed data", obsCount,
					      "(maybe specify weight='frequency' instead of numObs =",numObs,")"))
			}
		}
	}

	return(new("MxDataStatic", observed, means, type, as.numeric(numObs),
		observedStats, sort, primaryKey, weight, frequency, as.integer(verbose),
		as.logical(.parallel)))
}

setGeneric("preprocessDataForBackend", # DEPRECATED
	function(data, model, defVars, modeloptions) {
		return(standardGeneric("preprocessDataForBackend"))
	})

setGeneric("convertDataForBackend",
	function(data, model, flatModel) {
		return(standardGeneric("convertDataForBackend"))
	})

setGeneric("summarize",
	function(data) {
		return(standardGeneric("summarize"))
	})

setMethod("preprocessDataForBackend", signature("NonNullData"),
	  function(data, model, defVars, modeloptions) { data })

setMethod("convertDataForBackend", signature("NonNullData"),
	  function(data, model, flatModel) { data })

setMethod("convertDataForBackend", signature("MxDataStatic"),
	  function(data, model, flatModel) {
		  if (data@type == "cor") {
			  warning(paste("OpenMx does not yet correctly handle mxData(type='cor')",
					'standard errors and fit statistics.',
					'See Steiger (1980), "Tests for comparing elements of a correlation matrix".'))
		  }
		  if (data@type == "raw") {
			  if (is.matrix(data@observed) && is.integer(data@observed)) {
				  data@observed <- matrix(as.double(data@observed),
							  nrow=nrow(data@observed), ncol=ncol(data@observed),
							  dimnames=dimnames(data@observed))
			  }
			  if (is.data.frame(data@observed)) {
				  mapply(function(col, name) {
					  if (!is.factor(col)) return()
					  dups <- duplicated(levels(col))
					  if (any(dups)) {
						  stop(paste("Ordered factor column '",name,"' in model '",flatModel@name,
							     "' has more than 1 level with the same name: ",
							     omxQuotes(unique(levels(col)[dups])), sep=""))
					  }
				  }, data@observed, colnames(data@observed))
			  }
		  }
		  if (.hasSlot(data, 'primaryKey')) {
			  if (!is.na(data@primaryKey)) {
				  pk <- match(data@primaryKey, colnames(data@observed))
				  if (is.na(pk)) {
					  stop(paste("Primary key", omxQuotes(data@primaryKey),
						     "not found in", omxQuotes(data@name)))
				  }
				  data@primaryKey <- pk
			  }
		  }
		  if (.hasSlot(data, 'weight')) {
			  if (!is.na(data@weight)) {
				  wc <- match(data@weight, colnames(data@observed))
				  if (is.na(wc)) {
					  msg <- paste("Weight column", omxQuotes(data@weight),
						       "not found in observed data columns")
					  stop(msg, call.=FALSE)
				  }
				  data@weight <- wc
			  }
		  }
		  if (.hasSlot(data, 'frequency')) {
			  if (!is.na(data@frequency)) {
				  wc <- match(data@frequency, colnames(data@observed))
				  if (is.na(wc)) {
					  msg <- paste("Frequency column", omxQuotes(data@frequency),
						       "not found in observed data columns")
					  stop(msg, call.=FALSE)
				  }
				  data@frequency <- wc
			  }
		  }

		  data
	  })

setMethod("preprocessDataForBackend", signature("MxDataStatic"),
		function(data, model, defVars, modeloptions) {
			if(!is.null(data) && length(data@observedStats)) {
				obsStats <- data@observedStats
				if (!is.null(obsStats$thresholds)) {
					verifyThresholdNames(obsStats$thresholds, data@observed, model@name)
				}
			}
			data
		})

setMethod("summarize", signature("MxDataStatic"),
	  function(data) {
		  if (data@type != "raw") {
			  result <- list()
			  result[[data@type]] <- data@observed
			  if (!single.na(data@means)) {
				  result[['means']] <- data@means
			  }
			  return(result)
		  } else {
			  return(summary(data@observed))
		  }
	  })

setMethod("convertDataForBackend", signature("MxDataDynamic"),
	  function(data, model, flatModel) {
		  expNum <- match(data@expectation, names(flatModel@expectations))
		  if (any(is.na(expNum))) {
			  missing <- data@expectation[is.na(expNum)]
			  stop(paste("Cannot find expectation", omxQuotes(missing),
						"referenced by data in", model@name))
		  }
		  data@expectation <- expNum - 1L
		  data
	  })

setMethod("summarize", signature("MxDataDynamic"),
	  function(data) {
		  # Maybe we can do something smarter? TODO
		  c(dynamic=NA)
	  })

preprocessDatasets <- function(model, defVars, modeloptions) { # DEPRECATED
	if (!is.null(model@data)) {
		model@data <- preprocessDataForBackend(model@data, model, defVars, modeloptions)
	}
	if (length(model@submodels) > 0) {
		model@submodels <- lapply(model@submodels, preprocessDatasets,
					  defVars=defVars, modeloptions=modeloptions)
	}
	return(model)
}

convertDatasets <- function(datasets, model, flatModel) {
	lapply(datasets, convertDataForBackend, model=model, flatModel=flatModel)
}

checkNumericData <- function(data) {
	if(is.matrix(data@observed) && !is.double(data@observed)) {
		msg <- paste("The data object",
			omxQuotes(data@name), "contains an observed matrix that",
			"is not of type 'double'")
		stop(msg, call. = FALSE)
	}
}

getDataThresholdNames <- function(data) {
	if(data@type == "raw") {
		nCols <- ncol(data@observed)
		thresholdCols <- unlist(lapply(data@observed, is.ordered))
		if(sum(thresholdCols) > 0) {
			if(is.null(colnames(data@observed))) {
				stop("The data object", omxQuotes(data@name),
				"has ordinal columns, but does not have column names,",
				"which are required for threshold estimation.")
			}
			return(colnames(data@observed)[thresholdCols])
		}
	} else {
		obsStats <- data@observedStats
		if (is.null(obsStats$thresholds)) return(c())
		return(colnames(obsStats$thresholds))
	}
	return(c())
}

verifyCovarianceMatrix <- function(covMatrix, nameMatrix="observed") {
	if (is.null(covMatrix)) stop(paste("Covariance matrix",omxQuotes(nameMatrix),
		"is missing in action"))
	if(nrow(covMatrix) != ncol(covMatrix)) {
		msg <- paste("The", nameMatrix, "covariance matrix",
			"is not a square matrix - perhaps you meant type = 'raw' instead of 'cov'? ")
		stop(msg, call. = FALSE)
	}
	if (any(is.na(covMatrix))) {
		msg <- paste("The", nameMatrix, "covariance matrix",
			"contains NA values. Perhaps ensure you are excluding NAs in your cov() statement?")
		stop(msg, call. = FALSE)	
	}
	if (max(abs(covMatrix - t(covMatrix))) > 1e-9) {
		if(all.equal(covMatrix, t(covMatrix))){
			msg <- paste("The", nameMatrix, "covariance matrix",
				"is not a symmetric matrix,",
				" possibly due to rounding errors.\n",
				"Something like this would be appropriate:\n",
				"covMatrix[lower.tri(covMatrix)] = t(covMatrix)[lower.tri(t(covMatrix))]\n",
				"Where covMatrix is the name of your covariance data.",
				"Another option is \n round(covMatrix, 6)\n.")
		} else {
			msg <- paste("The", nameMatrix, "covariance matrix",
				"is not a symmetric matrix.\n",
				"Check what you are providing to mxData",
				"and perhaps try round(yourData, x) for x digits of precision.")
		}		
		stop(msg, call. = FALSE)
	}
	evalCov <- eigen(covMatrix)$values
	if (nameMatrix=="observed" & any( evalCov <= 0)) {
		msg <- paste("The", nameMatrix, "covariance matrix",
			"is not a positive-definite matrix:\n",
			"1 or more elements of eigen(covMatrix)$values ",
			"<= 0")
		stop(msg, call. = FALSE)
	}
	if (nameMatrix=="asymptotic" & any( evalCov < -1e-4)) {
		msg <- paste("The", nameMatrix, "covariance matrix",
			"is not a positive-semi-definite matrix:\n",
			"1 or more elements of eigen(covMatrix)$values ",
			"< 0")
		stop(msg, call. = FALSE)
	}
}

verifyCorrelationMatrix <- function(corMatrix) {
	if(nrow(corMatrix) != ncol(corMatrix)) {
		msg <- paste("The observed correlation matrix",
			"is not a square matrix")
		stop(msg, call. = FALSE)
	}
	if (any(is.na(corMatrix))) {
		msg <- paste("The observed correlation matrix",
			"contains NA values")
		stop(msg, call. = FALSE)	
	}
	if (!all(corMatrix == t(corMatrix))) {
		msg <- paste("The observed correlation matrix",
			"is not a symmetric matrix.",
			"Check what you are providing to mxData, and perhaps try using",
			"round(yourData, x) for x digits of precision.")
		stop(msg, call. = FALSE)
	}
	if (any(eigen(corMatrix)$values <= 0)) {
		msg <- paste("The observed correlation matrix",
			"is not a positive-definite matrix")
		stop(msg, call. = FALSE)
	}
	if (!all(diag(corMatrix) == 1)) {
		if (max(abs(diag(corMatrix) - 1)) < 1.0e-8) {
			msg <- paste("The observed correlation matrix",
				"has values that are near, but not equal to,",
				"1's along the diagonal")
		} else {
			msg <- paste("The observed correlation matrix",
				"is not 1's along the diagonal")
		}
		stop(msg, call. = FALSE)	
	}
}

generateDataList <- function(model) {
	retval <- lapply(model@submodels, generateDataList)
	names(retval) <- NULL
	retval <- unlist(retval)
	if (is.null(retval)) retval <- list()
	retval[[model@name]] <- model@data
	return(retval)
}

undoDataShare <- function(model, dataList) {
	model@data <- dataList[[model@name]]
	return(model)
}


displayMxData <- function(object) {
	cat("MxData", omxQuotes(object@name), '\n')
	cat("type :", omxQuotes(object@type), '\n')
	if (.hasSlot(object, 'primaryKey') && !is.na(object@primaryKey)) {
		cat("primary key :", omxQuotes(object@primaryKey), '\n')
	}
	cat("numObs :", omxQuotes(object@numObs), '\n')
	cat("observed : \n") 
	print(object@observed)
	if (length(object@means) == 1 && is.na(object@means)) {
		cat("means : NA \n")
	} else {
		cat("means : \n") 
		print(object@means)		
	}
	if (.hasSlot(object, 'observedStats') && length(object@observedStats)) {
		cat("observedStats : \n")
		print(object@observedStats)
	}
	invisible(object)
}

setMethod("print", "MxDataStatic", function(x,...) {
	displayMxData(x) 
})

setMethod("show", "MxDataStatic", function(object) {
	displayMxData(object) 
})

displayMxDataDynamic <- function(object) {
	cat("MxDataDynamic", omxQuotes(object@name), '\n')
	cat("type :", omxQuotes(object@type), '\n')
	cat("verbose :", omxQuotes(object@verbose), '\n')
	cat("expectation :", omxQuotes(object@expectation), '\n')
	invisible(object)
}

setMethod("print", "MxDataDynamic", function(x,...) {
	displayMxDataDynamic(x) 
})

setMethod("show", "MxDataDynamic", function(object) {
	displayMxDataDynamic(object) 
})
