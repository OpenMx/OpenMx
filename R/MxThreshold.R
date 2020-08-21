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

setClass(Class = "MxThreshold",
	representation = representation(
		variable = "character",
		nThresh  = "numeric",
		free     = "logical",
    values   = "numeric",
		labels   = "character",
		lbound   = "numeric",
		ubound   = "numeric"
))

setMethod("initialize", "MxThreshold",
	function(.Object, variable, nThresh, free, values, labels, lbound, ubound) {
		.Object@variable <- variable
		.Object@nThresh  <- nThresh
		.Object@free     <- free
		.Object@values   <- values
		.Object@labels   <- labels
		.Object@lbound   <- lbound
		.Object@ubound   <- ubound
		return(.Object)
	}
)

setMethod("$", "MxThreshold", imxExtractSlot)

setReplaceMethod("$", "MxThreshold",
	function(x, name, value) {
        stop("Changing threshold values directly is not recommended.  Please use the mxThreshold() function instead.")
	}
)

setMethod("names", "MxThreshold", slotNames)

##' mxNormalQuantiles
##'
##' Get quantiles from a normal distribution
##'
##' @param nBreaks the number of thresholds, or a vector of the number of thresholds
##' @param mean the mean of the underlying normal distribution
##' @param sd the standard deviation of the underlying normal distribution
##' @aliases omxNormalQuantiles
##' @return
##' a vector of quantiles
##' @examples
##' mxNormalQuantiles(3)
##' mxNormalQuantiles(3, mean=7)
##' mxNormalQuantiles(2, mean=1, sd=3)
mxNormalQuantiles <- function(nBreaks, mean=0, sd=1) {
	if(length(nBreaks) > 1) {
		return(unlist(mapply(mxNormalQuantiles, nBreaks=nBreaks, mean=mean, sd=sd)))
	}
	if(is.na(nBreaks)) {
		return(as.numeric(NA))
	}
	if(nBreaks < 0) {
		stop("Error in mxNormalQuantiles: You must request at least one quantile.")
	}
	return(qnorm(seq(0,1,length.out=nBreaks+2)[1:nBreaks+1], mean=mean, sd=sd))
}

omxNormalQuantiles <-mxNormalQuantiles

thresholdCheckLengths <- function(nThresh, starts, free, values, labels, lbounds, ubounds, varName) {

  if(single.na(starts)) {
    starts = rep(0, 5)
  }
  ends <- matrix(NA, 5, 2)

  ends[1,] <- thresholdCheckSingleLength(nThresh, starts[1], length(free), "free/fixed designations", varName)
  ends[2,] <- thresholdCheckSingleLength(nThresh, starts[2], length(values), "values", varName)
  ends[3,] <- thresholdCheckSingleLength(nThresh, starts[3], length(labels), "labels", varName)
  ends[4,] <- thresholdCheckSingleLength(nThresh, starts[4], length(lbounds), "lower bounds", varName)
  ends[5,] <- thresholdCheckSingleLength(nThresh, starts[5], length(ubounds), "upper bounds", varName)

  ends
}

thresholdCheckSingleLength <- function(nElements, startingPoint, len, lenName, varName) {

  if(nElements > len - startingPoint) {
		if(startingPoint != 0 && startingPoint %% len != 0) {
			stop(paste("This is tricky.  You've specified a set of",
				lenName, "for the mxThreshold function that does not",
				"quite line up. mxThreshold has consumed", startingPoint,
				"values from the list of", len, "but needs", nElements, "to complete",
				"the list of", lenName, "for", omxQuotes(varName)), call.=FALSE)
		}
    if(nElements %% len != 0) {
			stop(paste("mxThreshold() call will generate",
				nElements, "thresholds but you have specified",
				len, lenName, "for", omxQuotes(varName)), call. = FALSE)
		} else {
      startingPoint = 0
		}
	}
	return(c(startingPoint+1, startingPoint + nElements))
}

as.MxMatrix.MxThreshold <- function(threshold) {
  mxMatrix("Full", nrow=threshold@nThresh, ncol=1, free=threshold@free, values=threshold@values,
           labels=threshold@labels, lbound=threshold@lbound, ubound=threshold@ubound,
           dimnames=list(NULL, threshold@variable),condenseSlots=FALSE)
}

thresholdSubset <- function(var, a) {
  return(rep(var, length.out=a[2])[a[1]:a[2]])
}

splitThresholds <- function(variables, nThresh, free, values, labels, lbound, ubound) {

	nVars <- length(variables)
	nThresh <- rep(nThresh, length.out=nVars)
	nElements <- sum(nThresh)
	thresholds <- as.list(rep(NA, nVars))

  starts <- NA
  # This block is primarily to catch cases like nThresh=c(2, 4, 2), free=c(TRUE, TRUE, TRUE, FALSE)
  # That is, where the wrapping of an input vector doesn't mesh with the threshold count
  for(i in 1:nVars) {
    idx <- thresholdCheckLengths(nThresh[i], starts, free, values, labels, lbound, ubound, variables[i])
    thresholds[[i]] <- new("MxThreshold", variables[i], nThresh[i], thresholdSubset(free, idx[1,]),
                         as.numeric(thresholdSubset(values, idx[2,])),
                         as.character(thresholdSubset(labels, idx[3,])),
                         as.numeric(thresholdSubset(lbound, idx[4,])),
                         as.numeric(thresholdSubset(ubound, idx[5,])))
    starts <- idx[,2]
  }

  if(length(thresholds) == 1) {
    thresholds <- thresholds[[1]]
  }

	return(thresholds)
}

mxThreshold <- function(vars, nThresh=NA, free=FALSE, values=mxNormalQuantiles(nThresh),
                        labels=NA, lbound=NA, ubound=NA) {
  if(all.na(vars)) {
    stop("You must specify a variable name for which these thresholds should be applied.")
  }
	nVars = length(vars)
	if(single.na(nThresh)) {
		if(nVars > 0) {
			stop(paste("mxThreshold error: You must specify a list of numbers",
			" of levels if you specify more than one variable name."))
		} else if( all(is.na(values))) {
			stop(paste("mxThreshold error: You must specify either values or a",
			"number of levels for each variable in order to generate thresholds."))
		}
		nThresh = length(values)
	}
	repeats <- setdiff(vars, unique(vars))
	if(length(repeats) != 0) {
		stop(paste("mxThreshold error: You listed", omxQuotes(repeats),
		"more than once in your list of thresholded variables."))
	}

	nLevelsListed = length(nThresh)
	if(nVars < nLevelsListed) {
		stop(paste("mxThreshold error: You specified a number of levels for",
			 nLevelsListed, "variables, but only listed", nVars, "variables."))
	}
	if(!( nVars == nLevelsListed || nVars %% nLevelsListed==0)) {
			stop(paste("mxThreshold error: number of variables listed (",
			nVars, ") is not a multiple of the number of level counts (",
			nLevelsListed, ")."))
	}

	splitThresholds(vars, nThresh, free, values, labels, lbound, ubound)

}

mxMarginalProbit <- mxThreshold

generateDataThresholdColumns <- function(covarianceColumnNames, thresholdsMatrix) {
	covarianceLength <- length(covarianceColumnNames)
	thresholdColumns <- replicate(covarianceLength, as.integer(NA))
	thresholdLevels <- replicate(covarianceLength, as.integer(NA))

	thresholdColumnNames <- dimnames(thresholdsMatrix)[[2]]
	for(i in 1:length(thresholdColumnNames)) {
		oneThresholdName <- thresholdColumnNames[[i]]
		numThresholds <- sum(!is.na(thresholdsMatrix[,i]))
		columnIndex <- match(oneThresholdName, covarianceColumnNames)
		thresholdColumns[[columnIndex]] <- as.integer(i - 1)
		thresholdLevels[[columnIndex]] <- as.integer(numThresholds)
	}
	return(list(thresholdColumns, thresholdLevels))
}

verifyThresholds <- function(flatModel, model, labelsData, dataName, covNames, threshName) {
	if (single.na(threshName) || flatModel@unsafe) {
		return()
	}

	modelName <- flatModel@name
	object <- flatModel[[threshName]]
	if (is.null(object)) {
		stop(paste("The thresholds matrix/algebra", omxQuotes(threshName),
			"for model", omxQuotes(modelName),
			"does not exist."), call. = FALSE)
	} else if (is(object, "MxAlgebra") || is(object, "MxMatrix")) {
	} else {
		stop(paste("The thresholds entity", omxQuotes(threshName),
			"for model", omxQuotes(modelName),
			"is not a MxMatrix or MxAlgebra."), call. = FALSE)
	}

	tuple <- evaluateMxObject(threshName, flatModel, labelsData, new.env(parent = emptyenv()))
	thresholds <- tuple[[1]]
	observed <- flatModel@datasets[[dataName]]@observed
	dataType <- flatModel@datasets[[dataName]]@type

	threshNames <- verifyThresholdNames(thresholds, observed, modelName, observedThresholds=FALSE)

	dataThresh <- NA
	if(dataType == "acov") {
		dataThresh <- flatModel@datasets[[dataName]]@observedStats[['thresholds']]
	}

	for(i in 1:length(threshNames)) {
		tName <- threshNames[[i]]
		tColumn <- thresholds[,i]
		values = NA
		if(dataType=="raw") {
			observedColumn <- observed[,tName]
			if(!all(is.na(observedColumn))) {
				if (!is.ordered(observedColumn)) {
						stop(paste("In model",
							omxQuotes(modelName),
							"column",
							omxQuotes(tName),
							"is not an ordered factor.",
							"Use mxFactor() on this column."), call. = FALSE)
				}
				expectedThreshCount <- length(levels(observedColumn)) - 1
				if (nrow(thresholds) < expectedThreshCount) {
					stop(paste("In model",
						   omxQuotes(modelName),
						   "the number of thresholds in column",
						   omxQuotes(tName),
						   "is less than the (l - 1), where l is equal",
						   "to the number of levels in the ordinal",
						   "data. Use mxFactor() on this column."),
					     call. = FALSE)
				}
				values <- tColumn[1:expectedThreshCount]
				if (any(is.na(values))) {
					stop(paste("In model",
						omxQuotes(modelName),
						"I was expecting", expectedThreshCount,
						"thresholds in column", omxQuotes(tName),
						"of matrix/algebra",
						omxQuotes(simplifyName(threshName, modelName)),
						"but I hit NA values after only",
						length(values[!is.na(values)]), "thresholds.",
						"You need to increase the number of",
						"thresholds for", omxQuotes(tName),
						"and give them values other than NA"), call. = FALSE)
				}
			}
		} else if(dataType == "acov") {
			observedThresh <- dataThresh[,tName]
			if( !single.na(observedThresh)) {
				expectedThreshCount <- sum(!is.na(observedThresh))
				if (nrow(thresholds) < expectedThreshCount) {
					stop(paste("In model",
						omxQuotes(modelName),
						"the number of expected thresholds in column",
						omxQuotes(tName),
						"(", nrow(thresholds), "), is less than the number of thresholds",
						"observed in the data (", expectedThreshCount , "). If you use",
						"a prep function to set up your mxData object,",
						"be sure you ran mxFactor() on your data first.",
						"Otherwise, be sure all unused thresholds in your",
						"data thresholds argument are NA."),
						call. = FALSE)
				}
				values <- tColumn[1:expectedThreshCount]
				if (any(is.na(values))) {
					stop(paste("In model",
						omxQuotes(modelName),
						"I was expecting", expectedThreshCount,
						"thresholds in column", omxQuotes(tName),
						"of matrix/algebra",
						omxQuotes(simplifyName(threshName, modelName)),
						"but I hit NA values after only",
						length(values[!is.na(values)]), "thresholds.",
						"You need to increase the number of",
						"thresholds for", omxQuotes(tName),
						"and give them values other than NA"), call. = FALSE)
				}
			}
		}
		if(!single.na(values)) {
			sortValues <- sort(values, na.last = NA)
			if (!identical(sortValues, values)) {
				stop(paste("In model",
					omxQuotes(modelName),
					"the thresholds in column",
					omxQuotes(tName), "of matrix/algebra",
					omxQuotes(simplifyName(threshName, modelName)),
					"are not in ascending order.",
					"The current order is: ",
					omxQuotes(values), "and",
					"ascending order is: ",
					omxQuotes(sortValues),".",
					"Only the first", expectedThreshCount,
					"element(s) of this column are inspected."),
					call. = FALSE)
			}
		}
	}
}

verifyThresholdNames <- function(thresholds, observed, modelName=NA, observedThresholds=TRUE) {
	if(is.na(modelName)) {
		modelName = "[Model Name Unknown]"
	}

	if (is.null(dimnames(thresholds)) || is.null(dimnames(thresholds)[[2]])) {
		if(!observedThresholds) {
			stop(paste("The thresholds matrix/algebra",
			"for model", omxQuotes(modelName), "does not contain column names"),
			call. = FALSE)
		} else {
			stop(paste(
				"The observed thresholds matrix does not contain column names"),
				 call. = FALSE)
		}
	}
	if (is.null(dimnames(observed)) || is.null(dimnames(observed)[[2]])) {
		if(!observedThresholds) {
			stop(paste("The observed data frame for model",
				omxQuotes(modelName), "does not contain column names"),
				call. = FALSE)
		} else {
			stop(paste("The observed data does not contain column names"),
			call. = FALSE)
		}
	}
	threshNames <- dimnames(thresholds)[[2]]
	obsNames <- dimnames(observed)[[2]]
	missingNames <- setdiff(threshNames, obsNames)
	if (length(missingNames) > 0) {
		if(!observedThresholds) {
			stop(paste("The column name(s)", omxQuotes(missingNames),
				"appear in the thresholds matrix but not",
				"in the observed data.frame or matrix",
				"in model", omxQuotes(modelName)), call. = FALSE)
		} else {
			stop(paste("The column name(s)", omxQuotes(missingNames),
				"appear in the thresholds but not",
				"in the observed data."), call. = FALSE)
		}
	}
	return(threshNames)
}


displayThreshold <- function(object) {
	cat("mxThreshold", '\n')
	cat("$variable: ", omxQuotes(object@variable), '\n')
	cat("$nThresh", object@nThresh, '\n')
	cat("$values: ", object@values, '\n')
	cat("$free: ", object@free, '\n')
	cat("$labels: ", object@labels, '\n')
	cat("$lbound: ", object@lbound, '\n')
	cat("$ubound: ", object@ubound, '\n')
}

setMethod("print", "MxThreshold", function(x,...) { displayThreshold(x) })
setMethod("show", "MxThreshold", function(object) { displayThreshold(object) })
setAs("MxThreshold", "MxMatrix", function(from) { as.MxMatrix.MxThreshold(from)})

# --------------

#' An S4 base class for discrete marginal distributions
#'
#' @aliases DiscreteBase-class $,DiscreteBase-method $<-,DiscreteBase-method
#' @seealso \link{mxMarginalPoisson}, \link{mxMarginalNegativeBinomial}
DiscreteBase <- setClass(Class = "DiscreteBase",
	representation = representation(
		variable = "character",
		maxCount = "integer",
		free     = "logical",
		labels   = "character",
		lbound   = "numeric",
		ubound   = "numeric"))

setMethod("names", "DiscreteBase", slotNames)

setMethod("$", "DiscreteBase", imxExtractSlot)

setReplaceMethod("$", "DiscreteBase",
	function(x, name, value) {
    stop(paste("Changing",class(x), "values directly is not recommended."))
	}
)

setGeneric("getSpec", function(mxm) return(standardGeneric("getSpec")))

setClass(Class = "MxMarginalPoisson",
         contains = "DiscreteBase",
         representation = representation(
           zeroInf  = "numeric",
           lambda   = "numeric"))

setMethod("initialize", "MxMarginalPoisson",
          function(.Object, variable, maxCount, lambda, zeroInf, free,
                   labels, lbound, ubound) {
            .Object@variable <- variable
            .Object@maxCount <- as.integer(maxCount)
            .Object@lambda   <- lambda
            .Object@zeroInf  <- zeroInf
            .Object@free     <- as.logical(free)
            .Object@labels   <- as.character(labels)
            .Object@lbound   <- as.numeric(lbound)
            .Object@ubound   <- as.numeric(ubound)
            return(.Object)
          })

setMethod("getSpec", "MxMarginalPoisson", function(mxm) c(mxm@maxCount, 1))

setAs("MxMarginalPoisson", "MxMatrix", function(from) {
  mxMatrix("Full", nrow=2, ncol=1,
           free=from@free,
           values=c(from@zeroInf, from@lambda),
           labels=from@labels,
           lbound=from@lbound,
           ubound=from@ubound,
           dimnames=list(NULL, from@variable),
           condenseSlots=FALSE)
})

displayMarginalPoisson <- function(Ob) {
	cat("mxMarginalPoisson", '\n')
  for (sl in c("variable", "maxCount", "lambda", "zeroInf",
               "free", "labels", "lbound", "ubound")) {
    slname <- paste0("$", sl)
    cat(slname, ":", slot(Ob, sl), '\n')
  }
  invisible(Ob)
}

setMethod("print", "MxMarginalPoisson", function(x,...) { displayMarginalPoisson(x) })
setMethod("show", "MxMarginalPoisson", function(object) { displayMarginalPoisson(object) })

#' Indicator with marginal Poisson distribution
#'
#' @param vars character vector of manifest indicators
#' @param maxCount maximum observed count
#' @param lambda non-negative means
#' @param zeroInf zero inflation parameter in probability units
#' @param free logical vector indicating whether paremeters are free
#' @param labels character vector of parameter labels
#' @param lbound numeric vector of lower bounds
#' @param ubound numeric vector of upper bounds
#' @return a list of MxMarginPoisson obects
#' @aliases MxMarginalPoisson-class print,MxMarginalPoisson-method show,MxMarginalPoisson-method $,MxMarginalPoisson-method $<-,MxMarginalPoisson-method
mxMarginalPoisson <- function(vars, maxCount=NA, lambda, zeroInf=.01,
                      free=TRUE, labels=NA, lbound=0, ubound=c(1,NA))
{
  for (par in c('maxCount','lambda','zeroInf')) {
    if (length(get(par)) != length(vars) &&
          length(vars) %% length(get(par)) != 0) {
      stop(paste("Parameter",omxQuotes(par),"has wrong length"))
    }
  }
  for (par in c('free','labels','lbound','ubound')) {
    if (length(get(par)) > 1 && length(get(par))/2 != length(vars) &&
          length(vars) %% (length(get(par))/2) != 0) {
      stop(paste("Parameter",omxQuotes(par),"has wrong length"))
    }
  }

	poi <- vector('list', length(vars))
  vx <- 1
  for (xx in 1:length(vars)) {
    poi[[xx]] <-
      new("MxMarginalPoisson", vars[xx],
          maxCount[1 + (xx-1) %% length(maxCount)],
          lambda[1 + (xx-1) %% length(lambda)],
          zeroInf[1 + (xx-1) %% length(zeroInf)],
          sapply(0:1, function(x) free[1 + (vx+x-1) %% length(free)]),
          sapply(0:1, function(x) labels[1 + (vx+x-1) %% length(labels)]),
          sapply(0:1, function(x) lbound[1 + (vx+x-1) %% length(lbound)]),
          sapply(0:1, function(x) ubound[1 + (vx+x-1) %% length(ubound)]))
    vx <- vx + 2
  }
  if (length(poi) == 1) poi <- poi[[1]]
  poi
}

setClass(Class = "MxMarginalNegativeBinomial",
         contains = "DiscreteBase",
         representation = representation(
           zeroInf  = "numeric",
           size   = "numeric",
           prob = "MxOptionalNumeric",
           mu = "MxOptionalNumeric"))

setMethod("initialize", "MxMarginalNegativeBinomial",
          function(.Object, variable, maxCount, size, prob, mu, zeroInf, free,
                   labels, lbound, ubound) {
            .Object@variable <- variable
            .Object@maxCount <- as.integer(maxCount)
            .Object@size     <- size
            .Object@prob     <- prob
            .Object@mu       <- mu
            .Object@zeroInf  <- zeroInf
            .Object@free     <- as.logical(free)
            .Object@labels   <- as.character(labels)
            .Object@lbound   <- as.numeric(lbound)
            .Object@ubound   <- as.numeric(ubound)
            return(.Object)
          })

setMethod("getSpec", "MxMarginalNegativeBinomial",
          function(mxm) c(mxm@maxCount, ifelse(length(mxm@prob), 2, 3)))

setAs("MxMarginalNegativeBinomial", "MxMatrix", function(from) {
  if (length(from@mu) == 0) {
    v2 <- from@prob
  } else {
    v2 <- from@mu
  }
  mxMatrix("Full", nrow=3, ncol=1,
           free=from@free,
           values=c(from@zeroInf, from@size, v2),
           labels=from@labels,
           lbound=from@lbound,
           ubound=from@ubound,
           dimnames=list(NULL, from@variable),
           condenseSlots=FALSE)
})

displayMarginalNegativeBinomial <- function(Ob) {
	cat("mxMarginalNegativeBinomial", '\n')
  for (sl in c("variable", "maxCount", "size", "prob", "mu", "zeroInf",
               "free", "labels", "lbound", "ubound")) {
    slname <- paste0("$", sl)
    cat(slname, ":", slot(Ob, sl), '\n')
  }
  invisible(Ob)
}

setMethod("print", "MxMarginalNegativeBinomial", function(x,...) { displayMarginalNegativeBinomial(x) })
setMethod("show", "MxMarginalNegativeBinomial", function(object) { displayMarginalNegativeBinomial(object) })

#' Indicator with marginal Negative Binomial distribution
#'
#' @param vars character vector of manifest indicators
#' @param maxCount maximum observed count
#' @param size positive target number of successful trials
#' @param prob probability of success in each trial
#' @param mu alternative parametrization via mean
#' @param zeroInf zero inflation parameter in probability units
#' @param free logical vector indicating whether paremeters are free
#' @param labels character vector of parameter labels
#' @param lbound numeric vector of lower bounds
#' @param ubound numeric vector of upper bounds
#' @return a list of MxMarginPoisson obects
#' @aliases MxMarginalNegativeBinomial-class print,MxMarginalNegativeBinomial-method show,MxMarginalNegativeBinomial-method $,MxMarginalNegativeBinomial-method $<-,MxMarginalNegativeBinomial-method
mxMarginalNegativeBinomial <- function(vars, maxCount=NA, size, prob=c(), mu=c(), zeroInf=.01,
                      free=TRUE, labels=NA, lbound=NA, ubound=NA)
{
  if (!missing(prob) && !missing(mu)) stop("'prob' and 'mu' both specified")
  isMu <- !missing(mu)
  if (!isMu) {
    if (missing(lbound)) lbound <- c(0,0,0)
    if (missing(ubound)) ubound <- c(1,NA,1)
  } else {
    if (missing(lbound)) lbound <- c(0,0,NA)
    if (missing(ubound)) ubound <- c(1,NA,NA)
  }
  parList <- c('maxCount','size','zeroInf', ifelse(isMu,'mu','prob'))

  for (par in parList) {
    if (length(get(par)) != length(vars) &&
          length(vars) %% length(get(par)) != 0) {
      stop(paste("Parameter",omxQuotes(par),"has wrong length"))
    }
  }
  for (par in c('free','labels','lbound','ubound')) {
    if (length(get(par)) > 1 && length(get(par))/3 != length(vars) &&
          length(vars) %% (length(get(par))/3) != 0) {
      stop(paste("Parameter",omxQuotes(par),"has wrong length"))
    }
  }

	poi <- vector('list', length(vars))
  vx <- 1
  for (xx in 1:length(vars)) {
    poi[[xx]] <-
      new("MxMarginalNegativeBinomial", vars[xx],
          maxCount[1 + (xx-1) %% length(maxCount)],
          size[1 + (xx-1) %% length(size)],
          prob[1 + (xx-1) %% length(prob)],
          mu[1 + (xx-1) %% length(mu)],
          zeroInf[1 + (xx-1) %% length(zeroInf)],
          sapply(0:2, function(x) free[1 + (vx+x-1) %% length(free)]),
          sapply(0:2, function(x) labels[1 + (vx+x-1) %% length(labels)]),
          sapply(0:2, function(x) lbound[1 + (vx+x-1) %% length(lbound)]),
          sapply(0:2, function(x) ubound[1 + (vx+x-1) %% length(ubound)]))
    vx <- vx + 3
  }
  if (length(poi) == 1) poi <- poi[[1]]
  poi
}
