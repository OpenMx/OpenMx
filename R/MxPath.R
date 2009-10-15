#
#   Copyright 2007-2009 The OpenMx Project
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

# returns a list of paths
generatePath <- function(from, to, 
		all, arrows, values, free,
		labels, lbound, ubound) {
		if (single.na(to)) {
			to <- from
			loop <- TRUE
		} else {
			loop <- FALSE
		}
		from <- as.list(from)
		to <- as.list(to)
		if (all) {
			from <- rep(from, each=length(to))	
		}
		if(!is.null(labels)) { 
			lapply(labels, omxVerifyReference, 
				paste("path from", omxQuotes(from), "to", omxQuotes(to)))
		}
        pathCheckLengths(from, to, arrows, values, 
            free, labels, lbound, ubound, loop)
		result <- suppressWarnings(mapply(generateSinglePath, from, to, 
		arrows, values, free,
		labels, lbound, ubound, SIMPLIFY = FALSE))
		return(result)
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
	select <- (values != 0)
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

pathCheckVector <- function(value, valname, type) {
	if (!is.vector(value)) {
		stop(paste("The", omxQuotes(valname), 
			"argument to mxPath must be a",
			type, "vector."), call. = FALSE)
	}
}

mxPath <- function(from, to = NA, all = FALSE, arrows = 1, free = TRUE,
	values = NA, labels = NA, lbound = NA, ubound = NA) {
	if (all != TRUE && all != FALSE) {
		stop("The 'all' argument to mxPath must be either true or false.", 
			call. = FALSE)
	}
	pathCheckVector(from, 'from', 'character')
	pathCheckVector(to, 'to', 'character')
	pathCheckVector(arrows, 'arrows', 'numeric')
	pathCheckVector(free, 'free', 'logical')
	pathCheckVector(labels, 'labels', 'character')
	pathCheckVector(values, 'values', 'numeric')
	pathCheckVector(lbound, 'lbound', 'numeric')
	pathCheckVector(ubound, 'ubound', 'numeric')
	from <- as.character(from)
	to <- as.character(to)
	arrows <- as.numeric(arrows)
	free <- as.logical(free)
	values <- as.numeric(values)
	labels <- as.character(labels)
	lbound <- as.numeric(lbound)
	ubound <- as.numeric(ubound)
	if (single.na(values)) values <- NULL
	if (single.na(labels)) labels <- NULL
	if (single.na(lbound)) lbound <- NULL
	if (single.na(ubound)) ubound <- NULL
	generatePath(from, to, all, arrows, 
		values, free, labels, 
		lbound, ubound)
}
