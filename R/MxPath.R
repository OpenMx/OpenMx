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
	    connect = "character",
	    joinKey = "character"
))

setMethod("initialize", "MxPath",
	function(.Object, from, to, arrows, values,
		free, labels, lbound, ubound, connect, joinKey) {
		.Object@from <- from
		.Object@to <- to
		.Object@arrows <- arrows
		.Object@values <- values
		.Object@free <- free
		.Object@labels <- labels
		.Object@lbound <- lbound
		.Object@ubound <- ubound
		.Object@connect <- connect
		.Object@joinKey <- joinKey
		return(.Object)
	}
)

setMethod("$", "MxPath", imxExtractSlot)

setReplaceMethod("$", "MxPath",
	function(x, name, value) {
        stop("Changing path values directly is not recommended.  Please use the mxPath() function to set path values.")
	}
)

setMethod("names", "MxPath", slotNames)

expandPathConnect <- function(from, to, connect) {
	# now expand the paths to check for errors
	excludeBivariate <- FALSE
	excludeSelf      <- FALSE
	
	# interpret 'connect' argument
	if ((connect == "unique.pairs" ) || (connect == "unique.bivariate")){
		excludeBivariate <- TRUE
	}
	if ((connect == "all.bivariate") || (connect == "unique.bivariate")){
		excludeSelf <- TRUE
	}
	
	# if a variable is a connect = "single" then it does not need to be expanded
	if ((connect != "single")){ 
	
		from <- rep(from, each=length(to))
		to   <- rep(to, length(from)/length(to))

		exclude <- rep(FALSE, length(from))

		# if 'excluderedundant', then exclude b,a if a,b is present
		if (excludeBivariate){
			sortedPairs <- t(apply(matrix(c(from, to), ncol = 2), 1, sort))
			exclude <- exclude | duplicated(sortedPairs)
		}

		# if 'excludeself', then exclude x,x paths
		if (excludeSelf){
			exclude <- exclude | (from==to)
		}
		
		from <- from[!exclude]
		to   <- to[!exclude]		
	}
	return(list(from=from,to=to))
}

# returns a list of paths
generatePath <- function(from, to,
		connect, arrows, values, free,
		labels, lbound, ubound, joinKey) {

	# save exactly what the user typed to pass to mxModel for creation
	unalteredTo <- to
	unalteredFrom <- from

	# check if user employed the loop shortcut by only specifying from
	if (single.na(to)) {
		loop <- TRUE
		to <- from
	} else {
		loop <- FALSE
	}
	
	expanded <- expandPathConnect(from, to, connect)
	from <- expanded$from
	to   <- expanded$to

	# check for a missing to or from
	pathCheckToAndFrom(from, to)

	if(any(labels %in% "one")){
		warning("It is unwise to use the word 'one' as a label.\n",
		"That has a special meaning, as it used in 'from = \"one\", ' in means paths.\n",
		"See help(mxPath) to learn about mxPaths and labels")
	}
	# check the labels for illegal references
	lapply(labels, imxVerifyReference, -1)
	
	# check for length mismatches
	pathCheckLengths(from, to, arrows, values, free, labels, lbound, ubound, loop)
	
	# create a new MxPath in the MxModel
	return(new("MxPath", unalteredFrom, unalteredTo, arrows, values, free, labels, lbound, ubound, connect, joinKey))
}

pathCheckToAndFrom <- function(from, to){
	# check for a missing to or from
	if (any(is.na(from)) || any(is.na(to))) {
		stop("The \'from\' field or the \'to\' field contains an NA", call. = FALSE)
	}
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

##' imxIsPath
##'
##' This is an internal function exported for those people who know
##' what they are doing.
##'
##' @param value value
imxIsPath <- function(value) {
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
	if (!is.vector(value) || !check(value) || length(value) == 0) {
		stop(paste("The", omxQuotes(valname), 
			"argument to mxPath must be a",
			type, "vector of length > 0 in",
			deparse(width.cutoff = 400L, imxLocateFunction("mxPath"))
			), call. = FALSE)
	}
}

mxPath <- function(from, to = NA, 
	connect = c("single", "all.pairs", "unique.pairs", 
	            "all.bivariate", "unique.bivariate"), arrows = 1, 
	free = TRUE, values = NA, labels = NA, lbound = NA, ubound = NA, ..., joinKey=as.character(NA)) {
	if (missing(from)) {
		stop("The 'from' argument to mxPath must have a value.")
	}
	if (is.logical(connect)) {
		msg <- paste("The 'all' argument to mxPath ",
			"has been deprecated. It has been replaced ",
			"with the safer interface 'connect' in OpenMx 1.2. ",
			"See ?mxPath for more information.")
		# throw an error if 'all' has illegal value
		if ((length(connect) != 1) || single.na(connect)) {
			stop(msg)
		} else {
			warning(msg)
		}
	}
	garbageArguments <- list(...)
	allArgument <- garbageArguments[['all']]
	if (!is.null(allArgument)) {
		msg <- paste("The 'all' argument to mxPath ",
			"has been deprecated. It has been replaced ",
			"with the safer interface 'connect' in OpenMx 1.2. ",
			"See ?mxPath for more information.")
		# throw an error if 'all' has illegal value
		if (!(is.logical(allArgument) && 
			(length(allArgument) == 1) && 
			!single.na(allArgument))) {
			stop(msg)
		# throw an error if 'all' and 'connect' are both specified
		} else if (!identical(connect, 
			c("single", "all.pairs", "unique.pairs", 
			"all.bivariate", "unique.bivariate"))) {
			stop(msg)
		} else {
			warning(msg)
			connect <- allArgument
		}
	} else if (length(garbageArguments) > 0) {
		msg <- paste("mxPath does not accept values",
				"for the '...' argument.",
				"See ?mxPath for more information.")
   		stop(msg)
    }
	if (identical(connect, c("single", "all.pairs", "unique.pairs", 
	            "all.bivariate", "unique.bivariate"))) {
		connect <- "single"
	}
	# eliminate this test when "all" argument is eliminated
	if (is.logical(connect)) {
		if (connect) {
			connect <- "all.pairs"
		} else {
			connect <- "single"
		}
	} else {
		if (!(length(connect) == 1 && !single.na(connect) && 
			 is.character(connect) && (connect %in% 
				c("single", "all.pairs", "unique.pairs", 
	            "all.bivariate", "unique.bivariate")))) {
			msg <- paste("'connect' must be one of",
					"'single', 'all.pairs', 'unique.pairs',",
		            "'all.bivariate', or 'unique.bivariate'")
			stop(msg)
		}
		if (identical(connect, "all.pairs") && identical(arrows, 2)) {
			msg <- paste("'connect=all.pairs' argument cannot",
						"be used with 'arrows=2.',",
						"Please use 'connect=unique.pairs'.")
			stop(msg)
		}
		if (identical(connect, "all.bivariate") && identical(arrows, 2)) {
			msg <- paste("'connect=all.bivariate' argument cannot",
						"be used with 'arrows=2'.",
						"Please use 'connect=unique.bivariate'.")
			stop(msg)
		}
		if (!identical(connect, "single") && length(arrows) != 1) {
			msg <- paste("multiple values for the 'arrows' argument are",
					"not allowed when the 'connect' argument",
					"is not equal to 'single'")
			stop(msg)
		}
	}
	if (length(joinKey) > 1) {
		msg <- paste("cannot only joinKey a single foreign key, not",
			     length(joinKey))
		stop(msg)
	}
	if (!is.na(joinKey)) {
		if (any(arrows != 1)) {
			msg <- paste("between level join mappings can only use single",
				     "headed arrows")
			stop(msg)
		}
	}
	if (all.na(to)) { to <- as.character(to) }
	if (all.na(from)) { from <- as.character(from) }
	if (all.na(values)) { values <- as.numeric(values) }
	if (all.na(labels)) { labels <- as.character(labels) }
	if (all.na(lbound)) { lbound <- as.numeric(lbound) }
	if (all.na(ubound)) { ubound <- as.numeric(ubound) }
	if (all.na(connect)) { connect <- as.character(connect) } 	
	pathCheckVector(from, 'from', is.character, 'character')
	pathCheckVector(to, 'to', is.character, 'character')
	pathCheckVector(arrows, 'arrows', is.numeric, 'numeric')
	pathCheckVector(free, 'free', is.logical, 'logical')
	pathCheckVector(labels, 'labels', is.character, 'character')
	pathCheckVector(values, 'values', is.numeric, 'numeric')
	pathCheckVector(lbound, 'lbound', is.numeric, 'numeric')
	pathCheckVector(ubound, 'ubound', is.numeric, 'numeric')
	generatePath(from, to, connect, arrows,
		values, free, labels, 
		lbound, ubound, joinKey)
}

displayPath <- function(object) {
	cat("mxPath", '\n')
	cat("$from: ", omxQuotes(object@from), '\n')
	cat("$to: ", omxQuotes(object@to), '\n')
	cat("$arrows: ", object@arrows, '\n')
	cat("$values: ", object@values, '\n')
	cat("$free: ", object@free, '\n')
	cat("$labels: ", object@labels, '\n')
	cat("$lbound: ", object@lbound, '\n')
	cat("$ubound: ", object@ubound, '\n')
    cat("$connect: ", object@connect, '\n')
    cat("$joinKey: ", object@joinKey, '\n')
}

setMethod("print", "MxPath", function(x,...) { displayPath(x) })
setMethod("show", "MxPath", function(object) { displayPath(object) })
