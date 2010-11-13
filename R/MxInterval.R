#
#   Copyright 2007-2010 The OpenMx Project
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
	return(omxInterval(reference, lowerValue, upperValue))
}

omxInterval <- function(reference, lowerdelta, upperdelta) {
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

generateIntervalList <- function(flatModel, useIntervals, modelname, parameters, useNames = FALSE) {
	if (length(useIntervals) != 1 || 
		typeof(useIntervals) != "logical" || 
		is.na(useIntervals)) {
		stop(paste("'intervals' argument", 
			"must be TRUE or FALSE in",
			deparse(width.cutoff = 400L, sys.call(-1))), call. = FALSE)
	}
	if (!useIntervals) {
		return(list())
	}
	retval <- lapply(flatModel@intervals, generateIntervalListHelper, flatModel, modelname, parameters, useNames)
	names(retval) <- NULL
	retval <- unlist(retval, recursive = FALSE)
	return(retval)
}


makeIntervalReference <- function(entityNumber, row, col, lower, upper) {
	return(c(entityNumber, row - 1, col - 1, lower, upper))
}

generateIntervalListHelper <- function(interval, flatModel, modelname, parameters, useNames) {
	pnames <- names(parameters)
	reference <- interval@reference
	entity <- flatModel[[reference]]
	if(reference %in% pnames) {
		location <- parameters[[reference]][[3]]
		location[[1]] <- - location[[1]] - 1
		retval <- list()
		if (useNames) {
			retval[[reference]] <- c(reference, NA, NA, interval@lowerdelta, interval@upperdelta)
		} else {
			retval[[reference]] <- c(location, 
				interval@lowerdelta, interval@upperdelta)
		}
		return(retval)
	} else if (!is.null(entity)) {
		entityValue <- eval(substitute(mxEval(x, flatModel, compute=TRUE),
			list(x = as.symbol(reference))))
		rows <- nrow(entityValue)
		cols <- ncol(entityValue)
		if (is(entity, "MxMatrix")) {
			free <- entity@free
		} else {
			free <- matrix(TRUE, rows, cols)
		}
		entityNumber <- omxLocateIndex(flatModel, reference, 
			paste("confidence interval", reference))
		retval <- list()
		for(i in 1:rows) {
			for(j in 1:cols) {
				if (free[i, j]) {
					newName <- paste(reference, '[', i, ',', j, ']', sep = '')
					if (useNames) {
						retval[[newName]] <- c(reference, i, j, interval@lowerdelta, interval@upperdelta)
					} else {
						retval[[newName]] <- makeIntervalReference(
							entityNumber, i, j, interval@lowerdelta, interval@upperdelta)
					}
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
		if (is(entity, "MxMatrix")) {
			free <- entity@free
		} else {
			free <- matrix(TRUE, rows, cols)
		}
		entityNumber <- omxLocateIndex(flatModel, entityName, 
			paste("confidence interval", reference))
		retval <- list()
		for(i in rows) {
			for(j in cols) {
				if (free[i, j]) {
					newName <- paste(entityName, '[', i, ',', j, ']', sep = '')
					if (useNames) {
						retval[[newName]] <- c(entityName, i, j, interval@lowerdelta, interval@upperdelta)					
					} else {
						retval[[newName]] <- makeIntervalReference(
							entityNumber, i, j, interval@lowerdelta, interval@upperdelta)
					}
				}
			}
		}
		return(retval)
	} else {
		stop(paste("Unknown reference to", omxQuotes(reference),
			"detected in a confidence interval",
			"specification in model", omxQuotes(modelname), "in",
			deparse(width.cutoff = 400L, sys.call(-3))), call. = FALSE)
	}
}

displayInterval <- function(object) {
	cat("MxInterval", '\n')
	cat("@reference: ", object@reference, '\n')
	cat("@lowerdelta: ", object@lowerdelta, '\n')
	cat("@upperdelta: ", object@upperdelta, '\n')
}

setMethod("print", "MxInterval", function(x,...) { displayInterval(x) })
setMethod("show", "MxInterval", function(object) { displayInterval(object) })


omxParallelCI <- function(model, run = TRUE, suppressWarnings = TRUE, maxRetry = 5) {
	if(missing(model) || !is(model, "MxModel")) {
		stop("first argument must be a MxModel object")
	}
	if(length(model@output) == 0) {
		stop("'model' argument to omxParallelCI must be a fitted model")
	}
	namespace <- omxGenerateNamespace(model)
	flatModel <- omxFlattenModel(model, namespace)	
	parameters <- generateParameterList(flatModel)
	intervals <- generateIntervalList(flatModel, TRUE, model@name, parameters, useNames = TRUE)
	if (length(intervals) == 0) return(model)
	template <- model
	template@intervals <- list()
	modelname <- model@name
	container <- mxModel(paste(modelname, "container", sep = "_"))
	submodels <- list()
	for(i in 1:length(intervals)) {
		interval <- intervals[[i]]
		newmodel <- mxModel(paste("interval", i, sep = ""))		
		if (!is.na(interval[[4]])) {
			newmodel <- mxModel(newmodel, generateSubmodels(interval, template, modelname, "lower", i))
		}
		if (!is.na(interval[[5]])) {
			newmodel <- mxModel(newmodel, generateSubmodels(interval, template, modelname, "upper", i))
		}
		submodels <- c(submodels, newmodel)
	}
	container <- mxModel(container, submodels)
	if (!run) {
		return(container)
	}
	container <- mxRun(container, suppressWarnings = suppressWarnings)
	container <- rerunIntervals(container, NULL, 0, maxRetry, intervals, suppressWarnings)
	tableCI <- matrix(as.numeric(NA), length(intervals), 2)
	tableCodes <- matrix(0, length(intervals), 2)
	dimnames(tableCI) <- list(names(intervals), c('lbound', 'ubound'))
	dimnames(tableCodes) <- list(names(intervals), c('lbound', 'ubound'))
	for(i in 1:length(intervals)) {
		interval <- intervals[[i]]
		submodel <- container@submodels[[i]]
		lowerInterval <- findSubmodel(submodel, "lower")
		upperInterval <- findSubmodel(submodel, "upper")
		if (!is.null(lowerInterval)) {
			tableCI[i, 'lbound'] <- mxEval(ci, lowerInterval)[1,1]
			tableCodes[i, 'lbound'] <- lowerInterval@output$status[[1]]
		}
		if (!is.null(upperInterval)) {
			tableCI[i, 'ubound'] <- mxEval(ci, upperInterval)[1,1]
			tableCodes[i, 'ubound'] <- upperInterval@output$status[[1]]
		}
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

needsJitter <- function(model1, model2) {
	retval <- abs(model1@output$estimate - model2@output$estimate)
	return(all(retval <= 10^-4))
}

jitterMatrix <- function(matrix) {
	select <- matrix@free
	pvalues <- matrix@values[select]
	matrix@values[select] <- pvalues + 10^-4
	return(matrix)
}

jitterModel <- function(model) {
	model@matrices <- lapply(model@matrices, jitterMatrix)
	if(length(model@submodels) > 0) {
		model@submodels <- lapply(model@submodels, jitterModel)
	}
	return(model)
}

rerunIntervals <- function(container, previous, retry, maxRetry, intervals, suppressWarnings) {
	if (retry >= maxRetry) {
		return(container)
	}
	newmodel <- mxModel("container")
	for(i in 1:length(intervals)) {
                submodel <- container@submodels[[i]]
                lowerInterval <- findSubmodel(submodel, "lower")
                upperInterval <- findSubmodel(submodel, "upper")
                if (!is.null(lowerInterval) && 
		   lowerInterval@output$status[[1]] != 0) {
			if (is.null(previous)) {
				newmodel <- mxModel(newmodel, lowerInterval)
			} else {
				prevSubmodel <- previous@submodels[[i]]
				prevLower <- findSubmodel(prevSubmodel, "lower")
				if (needsJitter(lowerInterval, prevLower)) {
					newmodel <- mxModel(newmodel, jitterModel(lowerInterval))
				} else {
					newmodel <- mxModel(newmodel, lowerInterval)
				}
			}
                }
                if (!is.null(lowerInterval) && 
		   upperInterval@output$status[[1]] != 0) {
			if (is.null(previous)) {
				newmodel <- mxModel(newmodel, upperInterval)
			} else {
				prevSubmodel <- previous@submodels[[i]]
				prevUpper <- findSubmodel(prevSubmodel, "upper")
				if (needsJitter(upperInterval, prevUpper)) {
					newmodel <- mxModel(newmodel, jitterModel(upperInterval))
				} else {
					newmodel <- mxModel(newmodel, upperInterval)
				}
			}
		}
	}
	if (length(newmodel@submodels) == 0) return(container)
	newmodel <- mxRun(newmodel, suppressWarnings = suppressWarnings)
	nextcontainer <- container
	nextoutput <- nextcontainer@output
	nextoutput$frontendTime <- nextoutput$frontendTime + newmodel@output$frontendTime
	nextoutput$backendTime <- nextoutput$backendTime + newmodel@output$backendTime
	nextoutput$independentTime <- nextoutput$independentTime + newmodel@output$independentTime
	nextoutput$wallTime <- nextoutput$wallTime + newmodel@output$wallTime
	nextcontainer@output <- nextoutput
	for(i in 1:length(newmodel@submodels)) {
		submodel <- newmodel@submodels[[i]]
		nextcontainer[[submodel@name]] <- submodel
	}
	return(rerunIntervals(nextcontainer, container, retry + 1, maxRetry, intervals, suppressWarnings))
}

findSubmodel <- function(model, suffix) {
	submodels <- model@submodels
	subnames <- names(submodels)
	match <- grep(paste(suffix, "$", sep = ""), subnames)
	if (length(match) == 0) {
		return(NULL)
	} else {
		return(submodels[[match]])
	}
}


generateSubmodels <- function(interval, template, modelname, type, num) {
	name <- paste(modelname, "ci", num, type, sep="_")
	minimum <- mxEval(objective, template)
	template <- mxRename(template, paste(name, "child", sep = "_"))
	if (is.na(interval[[2]]) && is.na(interval[[3]])) {
		reference <- interval[[1]]
	} else {
		entity <- interval[[1]]
		components <- unlist(strsplit(entity, omxSeparatorChar, fixed = TRUE))			
		if (length(components) == 1) {
			entity <- paste(template@name, entity, sep = ".")
		} else {
			entity <- renameReference(entity, modelname, template@name)
		}
		reference <- paste("`", entity, "`", "[", interval[[2]], ",", interval[[3]], "]", sep = "")
	}
	ci <- eval(substitute(mxAlgebra(x, name = 'ci'), 
		list(x = parse(text = reference)[[1]])))
	if (type == "lower") {
		algebra <- eval(substitute(mxAlgebra((offset + minimum - value) ^ 2 + x, name = "algObjective"),
			list(x = parse(text = reference)[[1]],
				minimum = minimum[1,1],
				value = as.symbol(paste(template@name, "objective", sep = ".")), 
				offset = as.numeric(interval[[4]]))))
	} else if (type == "upper") {
		algebra <- eval(substitute(mxAlgebra((offset + minimum - value) ^ 2 - x, name = "algObjective"),
			list(x = parse(text = reference)[[1]],
				minimum = minimum[1,1],
				value = as.symbol(paste(template@name, "objective", sep = ".")), 
				offset = as.numeric(interval[[5]]))))	
	} else {
		stop(paste("Illegal type", type, "in generateSubmodels"))
	}
	objective <- mxAlgebraObjective("algObjective")
	model <- mxModel(name, algebra, objective, ci, template, independent = TRUE)
	return(model)
}
