#
#   Copyright 2007-2016 The OpenMx Project
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

##' omxNormalQuantiles
##'
##' Get quantiles from a normal distribution
##' 
##' @param nBreaks the number of thresholds, or a vector of the number of thresholds
##' @param mean the mean of the underlying normal distribution
##' @param sd the standard deviation of the underlying normal distribution
##' @return
##' a vector of quantiles
##' @examples
##' omxNormalQuantiles(3)
##' omxNormalQuantiles(3, mean=7)
##' omxNormalQuantiles(2, mean=1, sd=3)
omxNormalQuantiles <- function(nBreaks, mean=0, sd=1) {
	if(length(nBreaks) > 1) {
		return(unlist(mapply(omxNormalQuantiles, nBreaks=nBreaks, mean=mean, sd=sd)))
	}
	if(is.na(nBreaks)) {
		return(as.numeric(NA))
	}
	if(nBreaks < 0) {
		stop("Error in omxNormalQuantiles: You must request at least one quantile.")
	}
	return(qnorm(seq(0,1,length.out=nBreaks+2)[1:nBreaks+1], mean=mean, sd=sd))
}

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
    thresholds[i] <- new("MxThreshold", variables[i], nThresh[i], thresholdSubset(free, idx[1,]), 
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

mxThreshold <- function(vars, nThresh=NA, free=FALSE, values=NA, # normalQuantiles(nThresh)
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

generateThresholdColumns <- function(flatModel, model, labelsData, covarianceColumnNames, dataName, threshName) {
	covarianceLength <- length(covarianceColumnNames)
	thresholdColumns <- replicate(covarianceLength, as.integer(NA))
	thresholdLevels <- replicate(covarianceLength, as.integer(NA))
	if (single.na(threshName)) {
		return(list(thresholdColumns, thresholdLevels))
	}

	tuple <- evaluateMxObject(threshName, flatModel, labelsData, new.env(parent = emptyenv()))
	thresholds <- tuple[[1]]

	thresholdColumnNames <- dimnames(thresholds)[[2]]
	datasource <- flatModel@datasets[[dataName]]@observed
	for(i in 1:length(thresholdColumnNames)) {
		oneThresholdName <- thresholdColumnNames[[i]]
		if(flatModel@datasets[[dataName]]@type != 'acov'){
			numThresholds <- length(levels(datasource[, oneThresholdName])) - 1
		} else {
			numThresholds <- sum(!is.na(flatModel@datasets[[dataName]]@thresholds[, oneThresholdName]))
			# Note: this requires NA to be the indicator for an ingorable threshold
		}
		columnIndex <- match(oneThresholdName, covarianceColumnNames)
		thresholdColumns[[columnIndex]] <- as.integer(i - 1)
		thresholdLevels[[columnIndex]] <- as.integer(numThresholds)
	}
	return(list(thresholdColumns, thresholdLevels))
}

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
	if (single.na(threshName)) {
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
		dataThresh <- flatModel@datasets[[dataName]]@thresholds
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

factorize <- function(x, levels, labels, exclude, collapse) {
	x <- as.character(x)
	if (length(exclude) && all(!is.na(exclude))) {
		overlap <- match(exclude, levels)
		if (any(!is.na(overlap))) {
			msg <- paste("Factor levels and exclude vector are not disjoint; both contain",
				     omxQuotes(levels[overlap]))
			stop(msg)
		}
		x[which(!is.na(match(x, exclude)))] <- NA
	}
	noMatch <- !is.na(x) & is.na(match(x, levels))
	if (any(noMatch)) {
		msg <- paste("The following values are not mapped to factor levels and not excluded:",
			     omxQuotes(unique(x[noMatch])))
		stop(msg)
	}
	if (collapse) {
    corder <- order(labels)
    cLabels <- labels[corder]
    cLevels <- levels[corder]
	  dups <- duplicated(cLabels)
	  newLevels <- cLevels[!dups]
		notDup <- which(!dups)
		for (dx in which(dups)) {
			from <- cLevels[dx]
			to <- newLevels[findInterval(dx, notDup)]
			x[x==from] <- to
		}
    mask <- !duplicated(labels)
		levels <- levels[mask]
		labels <- labels[mask]
	} else {
	  dups <- duplicated(labels)
	  if (any(dups)) stop(paste("Duplicate labels and collapse=TRUE not specified:",
					  omxQuotes(unique(labels[dups]))))
	}

	f <- factor(x, levels, labels, exclude, ordered=TRUE)
	attr(f, 'mxFactor') <- TRUE
	f
}

mxFactor <- function(x = character(), levels, labels = levels, exclude = NA, ordered = TRUE, collapse=FALSE) {
	if(missing(levels)) {
		stop("the 'levels' argument is not optional")
	}
	if(!identical(ordered, TRUE)) {
		stop("the 'ordered' argument must be TRUE")
	}
	if (is.data.frame(x)) {
		if (is.list(levels)) {
			return(data.frame(mapply(factorize, x, levels, labels,
				MoreArgs=list(exclude = exclude, collapse=collapse), SIMPLIFY=FALSE),
				check.names = FALSE, row.names=rownames(x)))
		} else {
			return(data.frame(lapply(x, factorize, levels, labels, exclude, collapse),
				check.names = FALSE, row.names=rownames(x)))
		}
	} else if (is.matrix(x)) {
		stop(paste("argument 'x' to mxFactor()",
		"is of illegal type matrix,",
		"legal types are vectors or data.frames"))
	} else {
		return(factorize(x, levels, labels, exclude, collapse))
	}
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
