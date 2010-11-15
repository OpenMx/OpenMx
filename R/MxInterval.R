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

generateIntervalList <- function(flatModel, useIntervals, modelname, parameters) {
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
	retval <- lapply(flatModel@intervals, generateIntervalListHelper, flatModel, modelname, parameters)
	names(retval) <- NULL
	retval <- unlist(retval, recursive = FALSE)
	return(retval)
}


makeIntervalReference <- function(entityNumber, row, col, lower, upper) {
	return(c(entityNumber, row - 1, col - 1, lower, upper))
}



generateIntervalListHelper <- function(interval, flatModel, modelname, parameters) {
	pnames <- names(parameters)
	reference <- interval@reference
	entity <- flatModel[[reference]]
	if(reference %in% pnames) {
		location <- parameters[[reference]][[3]]
		location[[1]] <- - location[[1]] - 1
		retval <- list()
		retval[[reference]] <- c(location, 
			interval@lowerdelta, interval@upperdelta)
		return(retval)
	} else if (!is.null(entity)) {
		entityValue <- eval(substitute(mxEval(x, flatModel, compute=TRUE),
			list(x = as.symbol(reference))))
		rows <- nrow(entityValue)
		cols <- ncol(entityValue)
		if (is(entity, "MxMatrix")) {
			free <- entity@free
			if (omxSymmetricMatrix(entity)) {
				free[lower.tri(free, diag=FALSE)] <- FALSE
			}
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
					retval[[newName]] <- makeIntervalReference(
						entityNumber, i, j, interval@lowerdelta, interval@upperdelta)
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
					retval[[newName]] <- makeIntervalReference(
						entityNumber, i, j, interval@lowerdelta, interval@upperdelta)
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
