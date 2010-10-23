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

setClass(Class = "MxPath",
	representation = representation(
		from = "character",
		to = "character",
		arrows = "numeric",
		values = "numeric",
		free = "logical",
		labels = "character",
		lbound = "numeric",
		ubound = "numeric",
		excludeself = "logical"
))

setMethod("initialize", "MxPath",
	function(.Object, from, to, arrows, values,
		free, labels, lbound, ubound, excludeself) {
		.Object@from <- from
		.Object@to <- to
		.Object@arrows <- arrows
		.Object@values <- values
		.Object@free <- free
		.Object@labels <- labels
		.Object@lbound <- lbound
		.Object@ubound <- ubound
		.Object@excludeself <- excludeself
		return(.Object)
	}
)

# returns a list of paths
generatePath <- function(from, to, 
		all, arrows, values, free,
		labels, lbound, ubound, excludeself) {
	if (single.na(to)) {
		to <- from
		loop <- TRUE
	} else {
		loop <- FALSE
	}
	if (all) {
		from <- rep(from, each=length(to))
	}
	missingvalues <- is.na(values)
	values[missingvalues] <- 0
	if(!is.null(labels)) { 
		lapply(labels, omxVerifyReference, -1)
	}
	pathCheckLengths(from, to, arrows, values, 
		free, labels, lbound, ubound, loop)
	return(new("MxPath", from, to, arrows, values, free, labels, lbound, ubound, excludeself))
}

pathCheckLengths <- function(from, to, arrows, values, 
        free, labels, lbound, ubound, loop) {
    numPaths <- max(length(from), length(to))
    pathCheckSingleLength(numPaths, length(arrows), "arrows", from, to, loop)
    pathCheckSingleLength(numPaths, length(values), "values", from, to, loop)
    pathCheckSingleLength(numPaths, length(free), "free/fixed designations", from, to, loop)
    pathCheckSingleLength(numPaths, length(labels), "labels", from, to, loop)
    pathCheckSingleLength(numPaths, length(lbound), "lbounds", from, to, loop)
    pathCheckSingleLength(numPaths, length(ubound), "ubounds", from, to, loop)
}

pathCheckSingleLength <- function(numPaths, len, lenName, from, to, loop) {
    if (numPaths < len) {
    	if (loop) { to <- NA }
        stop(paste("mxPath() call will generate", 
            numPaths, "paths but you have specified",
            len, lenName, "with 'from' argument assigned to", omxQuotes(from),
            "and 'to' argument assigned to", omxQuotes(to)), call. = FALSE)
    }
}


generateSinglePath <- function(from, to, 
		arrows, values, free,
		labels, lbound, ubound) {
	result <- list()
	result[['from']] <- from
	result[['to']] <- to
	result[['arrows']] <- arrows[[1]]
	result[['values']] <- values[[1]]	
	result[['free']] <- free[[1]]
	result[['labels']] <- labels[[1]]	
	result[['lbound']] <- lbound[[1]]
	result[['ubound']] <- ubound[[1]]
	return(result)
}

omxIsPath <- function(value) {
	return(is.list(value) && 
		!is.null(value[['from']]) &&
		!is.null(value[['to']]))
}

matrixToPaths <- function(mxMatrix, arrows = c(1,2)) {
	values <- mxMatrix@values
	free <- mxMatrix@free
	labels <- mxMatrix@labels
	if (arrows == 2) {
		values[upper.tri(values)] <- 0
		free[upper.tri(free)] <- FALSE
		labels[upper.tri(labels)] <- as.character(NA)
	}
	select <- (values != 0) | (free) | (!is.na(labels))
	if (length(select) > 0) {
 	    rowFactors <- row(values, as.factor=TRUE)
	    colFactors <- col(values, as.factor=TRUE)	
		fromNames <- as.character(colFactors[select])
		toNames <- as.character(rowFactors[select])
		if (length(fromNames) > 0 && length(toNames) > 0) {
			return(mxPath(from = fromNames, to = toNames, arrows = arrows))
		}
	}
	return(list())
}

meansToPaths <- function(mxMatrix) {
	if (is.null(mxMatrix)) return(list())
	values <- mxMatrix@values
	free <- mxMatrix@free
	labels <- mxMatrix@labels
	select <- (values != 0) | (free) | (!is.na(labels))
	if (length(select) > 0) {
	    colFactors <- col(values, as.factor=TRUE)
		toNames <- as.character(colFactors[select])
		if (length(toNames) > 0) {
			return(mxPath(from = 'one', to = toNames, arrows = 1))
		}
	}
	return(list())
}

pathCheckVector <- function(value, valname, check, type) {
	if (!is.vector(value) || !check(value)) {
		stop(paste("The", omxQuotes(valname), 
			"argument to mxPath must be a",
			type, "vector."), call. = FALSE)
	}
}

mxPath <- function(from, to = NA, all = FALSE, arrows = 1, free = TRUE,
	values = NA, labels = NA, lbound = NA, ubound = NA, excludeself = FALSE) {
	if (missing(from)) {
		stop("The 'from' argument to mxPath must have a value.")
	}
	if (length(all) != 1 || !is.logical(all) || is.na(all)) {
		stop("The 'all' argument to mxPath must be either true or false.")
	}
	if (length(excludeself) != 1 || !is.logical(excludeself) || is.na(excludeself)) {
		stop("The 'excludeself' argument to mxPath must be either true or false.")
	}
	if (all.na(to)) { to <- as.character(to) }
	if (all.na(values)) { values <- as.numeric(values) }
	if (all.na(labels)) { labels <- as.character(labels) }
	if (all.na(lbound)) { lbound <- as.numeric(lbound) }
	if (all.na(ubound)) { ubound <- as.numeric(ubound) }	
	pathCheckVector(from, 'from', is.character, 'character')
	pathCheckVector(to, 'to', is.character, 'character')
	pathCheckVector(arrows, 'arrows', is.numeric, 'numeric')
	pathCheckVector(free, 'free', is.logical, 'logical')
	pathCheckVector(labels, 'labels', is.character, 'character')
	pathCheckVector(values, 'values', is.numeric, 'numeric')
	pathCheckVector(lbound, 'lbound', is.numeric, 'numeric')
	pathCheckVector(ubound, 'ubound', is.numeric, 'numeric')
	generatePath(from, to, all, arrows, 
		values, free, labels, 
		lbound, ubound, excludeself)
}

displayPath <- function(object) {
	cat("mxPath", '\n')
	cat("@from: ", omxQuotes(object@from), '\n')
	cat("@to: ", omxQuotes(object@to), '\n')
	cat("@arrows: ", object@arrows, '\n')
	cat("@values: ", object@values, '\n')
	cat("@free: ", object@free, '\n')
	cat("@labels: ", object@labels, '\n')
	cat("@lbound: ", object@lbound, '\n')
	cat("@ubound: ", object@ubound, '\n')
	cat("@excludeself: ", object@excludeself, '\n')
}

setMethod("print", "MxPath", function(x,...) { displayPath(x) })
setMethod("show", "MxPath", function(object) { displayPath(object) })
