#
#   Copyright 2007-2014 The OpenMx Project
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


setClassUnion("MxDataFrameOrMatrix", c("data.frame", "matrix"))

setClass(Class = "NonNullData")

setClass(Class = "MxDataStatic",
	 contains = "NonNullData",
	 representation = representation(
		observed = "MxDataFrameOrMatrix",
		means  = "matrix",
		type   = "character",
		numObs = "numeric",
		acov   = "matrix",
		thresholds = "matrix",
		thresholdColumns = "integer",
		thresholdLevels = "integer",
		indexVector = "integer",
		identicalDefVars = "integer",
		identicalMissingness = "integer",
		identicalRows = "integer",
		.isSorted = "logical",
		.needSort = "logical",
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
	function(.Object, observed, means, type, numObs, acov, thresholds, sort, name = "data") {
		.Object@observed <- observed
		.Object@means <- means
		.Object@type <- type
		.Object@numObs <- numObs
		.Object@acov <- acov
		.Object@thresholds <- thresholds
		.Object@name <- name
		.Object@.needSort <- sort
		.Object@.isSorted <- FALSE
		return(.Object)
	}
)

setMethod("initialize", "MxDataDynamic",
	function(.Object, type, expectation, name = "data") {
		.Object@type <- type
		.Object@expectation <- expectation
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
imxDataTypes <- c("raw", "cov", "cor", "sscp", "acov")

##' Create dynamic data
##'
##' @param type type of data
##' @param ...  Not used.  Forces remaining arguments to be specified by name.
##' @param expectation the name of the expectation to provide the data
##' @aliases
##' MxDataDynamic-class
mxDataDynamic <- function(type, ..., expectation) {
	garbageArguments <- list(...)
	if (length(garbageArguments) > 0) {
		stop("mxDataDynamic does not accept values for the '...' argument")
	}
	if (type != "cov") stop("Type must be set to 'cov'")
	return(new("MxDataDynamic", type, expectation))
}

mxData <- function(observed, type, means = NA, numObs = NA, acov=NA, thresholds=NA, ...,
		   sort=TRUE) {
	garbageArguments <- list(...)
	if (length(garbageArguments) > 0) {
		stop("mxData does not accept values for the '...' argument")
	}
	if (length(means) == 1 && is.na(means)) means <- as.numeric(NA)
	if (length(acov) == 1 && is.na(acov)) acov <- matrix(as.numeric(NA))
	if (length(thresholds) == 1 && is.na(thresholds)) thresholds <- matrix(as.numeric(NA))
	if (missing(observed) || !is(observed, "MxDataFrameOrMatrix")) {
		stop("Observed argument is neither a data frame nor a matrix")
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
	if (!is.vector(means) || !is.numeric(means)) {
		stop("Means argument must be of numeric vector type")
	}
	if (type != "raw" && is.na(numObs)) {
		stop("Number of observations must be specified for non-raw data, i.e., add numObs=XXX to mxData()")
	}
	if (type == "raw" && missing(numObs)) {
		numObs <- nrow(observed)
	}
	if (type == "cov") {
		verifyCovarianceMatrix(observed)
	}
	if (type == "cor") {
		verifyCorrelationMatrix(observed)
	}
	if (type == "acov") {
		verifyCovarianceMatrix(observed)
		verifyCovarianceMatrix(acov, nameMatrix="asymptotic")
		if ( !single.na(thresholds) ) {
			verifyThresholdNames(thresholds, observed)
		}
	}
	if (type != "acov") {
		if (any(!is.na(acov))) {
			stop("The acov argument to mxData is only allowed for data with type='acov'")
		}
		if (any(!is.na(thresholds))) {
			stop("The thresholds argument to mxData is only allowed for data with type='acov'")
		}
	}
	numObs <- as.numeric(numObs)
	lapply(dimnames(observed)[[2]], imxVerifyName, -1)
	meanNames <- names(means)
	means <- as.matrix(means)
	dim(means) <- c(1, length(means))
	colnames(means) <- meanNames
	return(new("MxDataStatic", observed, means, type, numObs, acov, thresholds, sort))
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

setMethod("preprocessDataForBackend", signature("MxDataStatic"),
	  function(data, model, defVars, modeloptions) {
		  if (data@.needSort) {
			  data <- sortRawData(data, defVars, model@name, modeloptions)
			  data <- convertIntegerColumns(data)
		  }

		  if(!is.null(data) && !single.na(data@thresholds)) {
			  verifyThresholdNames(data@thresholds, data@observed, model@name)
			  retval <- generateDataThresholdColumns(covarianceColumnNames=dimnames(data@observed)[[2]],
								 thresholdsMatrix=data@thresholds)
			  data@thresholdColumns <- retval[[1]]
			  data@thresholdLevels <- retval[[2]]
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
		  if (is.na(expNum)) stop(paste("Cannot find expectation", data@expectation,
						"referenced by data in", model@name))
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

sortRawData <- function(mxData, defVars, modelname, modeloptions) {
	if (is.null(mxData)) {
		return(mxData)
	}
	if (mxData@type != "raw") {
		return(mxData)	
	}
	if (!mxData@.needSort) stop("No sort needed")

	observed <- mxData@observed
	nosort <- as.character(modeloptions[['No Sort Data']])
	fullname <- paste(modelname, 'data', sep = '.')
	components <- unlist(strsplit(fullname, imxSeparatorChar, fixed = TRUE))
	modelname <- components[[1]]
	if ((length(observed) == 0) || (modelname %in% nosort)) {
		mxData@indexVector <- as.integer(NA)
		mxData@identicalDefVars <- as.integer(NA)
		mxData@identicalMissingness <- as.integer(NA)
		mxData@identicalRows <- as.integer(NA)
	} else {
		observedNames <- colnames(observed)	
		if (is.null(observedNames)) {
			msg <- paste("The raw data set in model", omxQuotes(modelname),
				"does not contain column names")
			stop(msg, call. = FALSE)
		}
		if (length(defVars) > 0) {
			defkeys <- names(imxFilterDefinitionVariables(defVars, fullname))
			defkeys <- sapply(defkeys, function(x) {
				unlist(strsplit(x, imxSeparatorChar, fixed = TRUE))[[3]]
			})
			names(defkeys) <- NULL
			defkeys <- defkeys[defkeys %in% observedNames]
			if (length(defkeys) == 0) {
				defkeys <- character()
			}
			defindex <- match(defkeys, observedNames)
			otherindex <- setdiff(1:length(observedNames), defindex)
		} else {
			defkeys <- character()
			otherindex <- c(1:length(observedNames))
		}
		nacount <- sapply(otherindex, function(x) { sum(is.na(observed[,x])) })
		otherindex <- otherindex[order(nacount, decreasing=TRUE)]
		defvectors <- lapply(defkeys, function(x) {observed[,x] })
		othervectorsNA <- lapply(otherindex, function(x) {!is.na(observed[,x]) })
		othervectors <- lapply(otherindex, function(x) {observed[,x] })
		args <- c(defvectors, othervectorsNA, othervectors, 'na.last'=FALSE)
		indexVector <- do.call('order', args)
		sortdata <- observed[indexVector,,drop=FALSE]
		mxData@observed <- sortdata
		selectMissing <- is.na(sortdata)
		selectDefvars <- sortdata[, defkeys, drop=FALSE]
		threeVectors <- .Call(findIdenticalRowsData, sortdata,
			selectMissing, selectDefvars, sum(nacount)==0,
			length(selectDefvars) == 0, PACKAGE = "OpenMx")
		mxData@indexVector <- indexVector - 1L
		mxData@identicalRows <- threeVectors[[1]]
		mxData@identicalMissingness <- threeVectors[[2]]
		mxData@identicalDefVars <- threeVectors[[3]]
	}
	mxData@.needSort <- FALSE
	mxData@.isSorted <- TRUE
	return(mxData)
}

convertIntegerColumns <- function(mxData) {
	if (is.null(mxData)) return(mxData)
	if (is.data.frame(mxData@observed)) {
		dimnames <- dimnames(mxData@observed)
		mxData@observed <- data.frame(lapply(mxData@observed, convertIntegerSingleColumn))
		dimnames(mxData@observed) <- dimnames
	}
	mxData@numObs <- as.numeric(mxData@numObs)
	return(mxData)
}

convertIntegerSingleColumn <- function(column) {
	if(all(is.na(column)) || (!is.factor(column) && is.integer(column))) {
		return(as.double(column))
	} else {
		return(column)
	}
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
		thresholdCols <- rep(FALSE, nCols)
		for(i in 1:nCols) {
			if(!is.null(attr(data@observed[,i], "mxFactor"))) {
				thresholdCols[i] <- TRUE
			}
		}
		if(sum(thresholdCols) > 0) {
			if(is.null(colnames(data@observed))) {
				stop("The data object", omxQuotes(data@name),
				"has mxFactor columns, but does not have column names,",
				"which are required for threshold estimation.")
			}
			return(colnames(data@observed)[thresholdCols])
		}
	} else if(data@type == "acov") {
		return(dimnames(data@thresholds)[[2]])
	}
	return(c())
}

# MCN - Skip check and leave it to SADMVN
  checkNumberOrdinalColumns <- function(data) {
  	if(is.data.frame(data@observed)) {
  		count <- sum(sapply(data@observed, is.factor))
  		if (count > 20) {
#  			msg <- paste("The data object",
#  				omxQuotes(data@name), "contains", count,
#  				"ordered factors but our ordinal integration",
#  				"implementation has a limit of 20 ordered factors.")
#  			stop(msg, call. = FALSE)
  		}
  	}
  }
  
verifyCovarianceMatrix <- function(covMatrix, nameMatrix="observed") {
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
	if (!all(covMatrix == t(covMatrix))) {
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
	cat("numObs :", omxQuotes(object@numObs), '\n')
	cat("observed : \n") 
	print(object@observed)
	if (length(object@means) == 1 && is.na(object@means)) {
		cat("means : NA \n")
	} else {
		cat("means : \n") 
		print(object@means)		
	}
	if (length(object@acov) == 1 && is.na(object@acov)) {
		cat("acov : NA \n")
	} else {
		cat("acov : \n")
		print(object@acov)
	}
	if (length(object@thresholds) == 1 && is.na(object@thresholds)) {
		cat("thresholds : NA \n")
	} else {
		cat("thresholds : \n")
		print(object@thresholds)
	}
	invisible(object)
}

setMethod("print", "MxDataStatic", function(x,...) {
	displayMxData(x) 
})

setMethod("show", "MxDataStatic", function(object) {
	displayMxData(object) 
})
