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

##' MxInterval
##'
##' @description
##' This is an internal class and should not be used directly.
##'
##' @aliases
##' $,MxInterval-method
##' $<-,MxInterval-method
##' print,MxInterval-method
##' show,MxInterval-method
##' @seealso
##' \link{mxCI}
##' @rdname MxInterval-class
setClass(Class = "MxInterval",
	representation = representation(
		reference = "character",
		lowerdelta = "numeric",
		upperdelta = "numeric"
	))

setMethod("initialize", "MxInterval",
	function(.Object, reference, lowerdelta, upperdelta) {
		.Object@reference <- reference
		.Object@lowerdelta <- lowerdelta
		.Object@upperdelta <- upperdelta
		return(.Object)
	}
)


createNewInterval <- function(reference, lowerdelta, upperdelta) {
	return(new("MxInterval", reference, lowerdelta, upperdelta))
}

setMethod("$", "MxInterval", imxExtractSlot)

setReplaceMethod("$", "MxInterval",
	function(x, name, value) {
		return(imxReplaceSlot(x, name, value, check=TRUE))
	}
)

setMethod("names", "MxInterval", slotNames)

expandSingleInterval <- function(interval) {
	references <- interval@reference
	if (length(references) == 1) {
		return(interval)
	} else {
		return(lapply(references, createNewInterval, 
			interval@lowerdelta, interval@upperdelta))
	}
}

expandIntervals <- function(intervals) {
	if (length(intervals) == 0) {
		return(intervals)
	}
	retval <- lapply(intervals, expandSingleInterval)
	retval <- unlist(retval)
	return(retval)
}

mxCI <- function(reference, interval = 0.95, type = c('both', 'lower', 'upper')) {
	if (!is.numeric(interval) || interval < 0 || interval > 1) {
		stop("'interval' must be a numeric value between 0 and 1")
	}
	if (identical(type, c('both', 'lower', 'upper'))) {
		type <- 'both'
	}
	if (!is.character(type) || length(type) != 1 || !(type %in% c('both', 'lower', 'upper'))) {
		stop("'type' must be either 'both' or 'lower' or 'upper'")
	}
	if (type == 'both') {
		lowerValue <- qchisq(interval, 1)
		upperValue <- qchisq(interval, 1)
	} else if (type == 'lower') {
		lowerValue <- qchisq(interval, 1)
		upperValue <- as.numeric(NA)
	} else if (type == 'upper') {
		lowerValue <- as.numeric(NA)
		upperValue <- qchisq(interval, 1)
	}
	return(confidenceIntervalHelper(reference, lowerValue, upperValue))
}

confidenceIntervalHelper <- function(reference, lowerdelta, upperdelta) {
	if (single.na(lowerdelta)) { lowerdelta <- as.numeric(NA) }
	if (single.na(upperdelta)) { upperdelta <- as.numeric(NA) }
	if (!is.character(reference) || length(reference) < 1 || any(is.na(reference))) {
		stop("'reference' argument must be a character vector")
	}
	if (!is.na(lowerdelta) && (!is.numeric(lowerdelta) || length(lowerdelta) != 1 || lowerdelta < 0)) {
		stop("'lowerdelta' argument must be a non-negative numeric value")
	}
	if (!is.na(upperdelta) && (!is.numeric(upperdelta) || length(upperdelta) != 1 || upperdelta < 0)) {
		stop("'upperdelta' argument must be a non-negative numeric value")
	}
	retval <- createNewInterval(reference, lowerdelta, upperdelta)
	return(retval)
}

modelAddIntervals <- function(model, intervals) {
	if (length(intervals) == 0) {
		return(model)
	}
	iNames <- names(intervals)
	for(i in 1:length(intervals)) {
		model@intervals[[iNames[[i]]]] <- intervals[[i]]
	}
	return(model)
}

modelRemoveIntervals <- function(model, intervals) {
	if (length(intervals) == 0) {
		return(model)
	}
	iNames <- names(intervals)
	for(i in 1:length(intervals)) {		
		model@intervals[[iNames[[i]]]] <- NULL
	}
	return(model)
}

generateIntervalList <- function(flatModel, modelname, parameters, labelsData) {
	retval <- lapply(flatModel@intervals, generateIntervalListHelper, 
		flatModel, modelname, parameters, labelsData)
	names(retval) <- NULL
	retval <- unlist(retval, recursive = FALSE)
	return(retval)
}


makeIntervalReference <- function(entityNumber, row, col, lower, upper) {
	return(c(entityNumber, row - 1, col - 1, lower, upper))
}

generateIntervalListHelper <- function(interval, flatModel, modelname, 
			parameters, labelsData) {

	reference <- interval@reference
					# length(reference)==1 because of expandSingleInterval
	pindex <- match(reference, names(parameters))
	if (!is.na(pindex)) {
		location <- parameters[[pindex]][[5]]
		location[[1]] <- - location[[1]] - 1
		retval <- list()
		retval[[reference]] <-
		    c(location, interval@lowerdelta, interval@upperdelta)
		return(retval)
	}
	entity <- flatModel[[reference]]
	retval <- list()
	if (!is.null(entity)) {
		entityNumber <- imxLocateIndex(flatModel, reference, 
					       paste("confidence interval", reference))
		if (is(entity, "MxAlgebra")) {
			newName <- reference
			retval[[newName]] <- makeIntervalReference(entityNumber, 0L, 0L, 
								   interval@lowerdelta, interval@upperdelta)
		} else {
		  if(!any(entity@free)){return(retval)}
			rows <- nrow(entity)
			cols <- ncol(entity)
			free <- entity@free
			for(i in 1:rows) {
				for(j in 1:cols) {
					if (!free[i, j]) next
					if(!is.na(entity@labels[i,j]) ){
						newName <- entity@labels[i,j]
					} else {
						newName <- paste(reference, '[', i, ',', j, ']', sep = '')
					}
					retval[[newName]] <- makeIntervalReference(entityNumber, i, j, 
										   interval@lowerdelta, interval@upperdelta)
				}
			}
		}
		return(retval)
	} else if (hasSquareBrackets(reference)) {
		components <- splitSubstitution(reference)
		entityName <- components[[1]]
		rows <- eval(parse(text = components[[2]]), envir = globalenv())
		cols <- eval(parse(text = components[[3]]), envir = globalenv())
		entityValue <- eval(substitute(mxEval(x, flatModel, compute=TRUE),
			list(x = as.symbol(entityName))))
		if (is.null(rows)) {
			rows <- 1:nrow(entityValue)
		}
		if (is.null(cols)) {
			cols <- 1:ncol(entityValue)
		}
		entity <- flatModel[[entityName]]
		entityNumber <- imxLocateIndex(flatModel, entityName, 
			paste("confidence interval", reference))
		retval <- list()
		for(i in rows) {
			for(j in cols) {
				newName <- paste(entityName, '[', i, ',', j, ']', sep = '')
				retval[[newName]] <- makeIntervalReference(entityNumber, i, j, 
									   interval@lowerdelta, interval@upperdelta)
			}
		}
		return(retval)
	} else {
		for (entityName in names(flatModel@matrices)) {
			entity <- flatModel[[entityName]]
			free <- entity@free
			rows <- nrow(entity)
			cols <- ncol(entity)
			for(i in 1:rows) {
				for(j in 1:cols) {
					if (free[i, j]) next
					label <- entity@labels[i,j]
					if (is.na(label) || label != reference) next
					entityNumber <- imxLocateIndex(flatModel, entityName, 
								       paste("confidence interval", reference))
					retval[[label]] <- makeIntervalReference(entityNumber, i, j, 
										 interval@lowerdelta, interval@upperdelta)
					return(retval)
				}
			}
		}
		stop(paste("Unknown reference to", omxQuotes(reference),
			"detected in a confidence interval",
			"specification in model", omxQuotes(modelname), 
			"\nYou should check spelling (case-sensitive), and also addressing the right model: to refer to an algebra", 
			"\nSee help(mxCI) to see how to refer to an algebra in a submodel.\n", 
			"FYI, I got as far as: ", deparse(width.cutoff = 400L, sys.call(-3))), call. = FALSE)
	}
}

displayInterval <- function(object) {
	cat("MxInterval", '\n')
	cat("$reference: ", object@reference, '\n')
	cat("$lowerdelta: ", object@lowerdelta, '\n')
	cat("$upperdelta: ", object@upperdelta, '\n')
}

setMethod("print", "MxInterval", function(x,...) { displayInterval(x) })
setMethod("show", "MxInterval", function(object) { displayInterval(object) })

expandConfidenceIntervals <- function(model, intervals) {
	expansion <- lapply(intervals, expandConfidenceIntervalsHelper, model)
	return(unlist(expansion))
}

expandConfidenceIntervalsHelper <- function(interval, model) {
	reference <- interval@reference
	lowerdelta <- interval@lowerdelta
	upperdelta <- interval@upperdelta
	entity <- model[[reference]]
	retval <- list()
	if(!is.null(entity)) {
		entityValue <- eval(substitute(mxEval(x, model, compute=TRUE),
			list(x = as.symbol(reference))))
		rows <- nrow(entityValue)
		cols <- ncol(entityValue)
		if (is(entity, "MxMatrix")) {
			free <- entity@free
			if (imxSymmetricMatrix(entity)) {
				free[lower.tri(free, diag=FALSE)] <- FALSE
			}
		} else {
			free <- matrix(TRUE, rows, cols)
		}
		for(i in 1:rows) {
			for(j in 1:cols) {
				if (free[i, j]) {
					newName <- paste(reference, '[', i, ',', j, ']', sep = '')
					retval <- c(retval, new("MxInterval", newName, lowerdelta, upperdelta))
				}
			}
		}
	} else if (hasSquareBrackets(reference)) {
		components <- splitSubstitution(reference)
		entityName <- components[[1]]
		rows <- eval(parse(text = components[[2]]), envir = globalenv())
		cols <- eval(parse(text = components[[3]]), envir = globalenv())
		entityValue <- eval(substitute(mxEval(x, model, compute=TRUE),
			list(x = as.symbol(entityName))))
		if (is.null(rows)) {
			rows <- 1:nrow(entityValue)
		}
		if (is.null(cols)) {
			cols <- 1:ncol(entityValue)
		}
		entity <- model[[entityName]]
		if (is.null(rows)) {
			rows <- 1:nrow(entityValue)
		}
		if (is.null(cols)) {
			cols <- 1:ncol(entityValue)
		}
		entity <- model[[entityName]]
		if (is(entity, "MxMatrix")) {
			free <- entity@free
			if (imxSymmetricMatrix(entity)) {
				free[lower.tri(free, diag=FALSE)] <- FALSE
			}
		} else {
			free <- matrix(TRUE, rows, cols)
		}
		for(i in rows) {
			for(j in cols) {
				if (free[i, j]) {
					newName <- paste(entityName, '[', i, ',', j, ']', sep = '')
					retval <- c(retval, new("MxInterval", newName, lowerdelta, upperdelta))
				}
			}
		}
	} else {
		retval <- interval
	}
	return(retval)
}

removeAllIntervals <- function(model) {
	model@intervals <- list()
	model@submodels <- lapply(model@submodels, removeAllIntervals)
	return(model)
}

##' omxParallelCI
##'
##' Create parallel models for parallel confidence intervals
##' 
##' @param model an MxModel with confidence intervals in it
##' @param run whether to run the model or just return the parallelized interval models
##' @return
##' an MxModel object
##' @examples
##' require(OpenMx)
##' data(demoOneFactor)
##' manifests <- names(demoOneFactor)
##' latents <- c("G")
##' factorModel <- mxModel("One Factor", 
##'                       type="RAM",
##'                       manifestVars=manifests, 
##'                       latentVars=latents,
##'                       mxPath(from=latents, to=manifests),
##'                       mxPath(from=manifests, arrows=2),
##'                       mxPath(from=latents, arrows=2, free=FALSE, values=1.0),
##'                       mxData(observed=cov(demoOneFactor), type="cov",
##'                       numObs=500),
##'      # add confidence intervals for free params in A and S matrices
##'                       mxCI(c('A', 'S')))
##' factorRun <- mxRun(factorModel)
##' factorCI <- omxParallelCI(factorRun) # Run CIs in parallel
omxParallelCI <- function(model, run = TRUE) {
	if(missing(model) || !is(model, "MxModel")) {
		stop("first argument must be a MxModel object")
	}
	if(length(model@output) == 0) {
		stop("'model' argument to omxParallelCI must be a fitted model")
	}
	namespace <- imxGenerateNamespace(model)
	flatModel <- imxFlattenModel(model, namespace)
	intervals <- flatModel@intervals
	if (length(intervals) == 0) return(model)
	intervals <- expandConfidenceIntervals(model, intervals)
	template <- model
	template <- removeAllIntervals(template)
	modelname <- model@name
	container <- mxModel(paste(modelname, "container", sep = "_"))
	submodels <- list()
	for(i in 1:length(intervals)) {
		interval <- intervals[[i]]
		newmodel <- mxModel(template, interval, independent = TRUE)
		newmodel <- mxRename(newmodel, paste("interval", i, sep = ""))
		newmodel <- mxOption(newmodel, "Number of Threads", 1)
		submodels <- c(submodels, newmodel)
	}
	container <- mxModel(container, submodels)
	if (!run) {
		return(container)
	}
	container <- mxRun(container, intervals = TRUE, suppressWarnings = TRUE)
	tableCI <- matrix(as.numeric(NA), 0, 3)
	tableCodes <- matrix(0, 0, 2)
	dimnames(tableCI) <- list(NULL, c('lbound', 'estimate', 'ubound'))
	dimnames(tableCodes) <- list(NULL, c('lbound', 'ubound'))
	submodels <- container@submodels
	for(i in 1:length(submodels)) {
		submodel <- submodels[[i]]
		submodel <- mxRename(submodel, modelname)
		tableCI <- rbind(tableCI, submodel@output$confidenceIntervals)
		tableCodes <- rbind(tableCodes, submodel@output$confidenceIntervalCodes)
	}
	model@output$confidenceIntervals <- tableCI
	model@output$confidenceIntervalCodes <- tableCodes
	model@output$frontendTime <- container@output$frontendTime
	model@output$backendTime <- container@output$backendTime
	model@output$independentTime <- container@output$independentTime
	model@output$wallTime <- container@output$wallTime
	model@output$timestamp <- container@output$timestamp
	model@output$cpuTime <- container@output$cpuTime
	return(model)
}
