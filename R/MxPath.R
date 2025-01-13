#
#   Copyright 2007-2021 by the individuals mentioned in the source code history
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
    joinKey = "character",
    step="MxOptionalInteger",
		strictUnigraph="logical"
))

setMethod("initialize", "MxPath",
	function(.Object, from, to, arrows, values,
		free, labels, lbound, ubound, connect, joinKey, step, strictUnigraph) {
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
    .Object@step <- step
    .Object@strictUnigraph <- strictUnigraph
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
		labels, lbound, ubound, joinKey, step, strictUnigraph) {

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
	if (max(length(from), length(to)) == 0 && length(values) <= 1 &&
		    length(free) <= 1 && length(labels) <= 1 && length(lbound) <= 1 &&
		    length(ubound) <= 1 && length(joinKey) <= 1 && length(strictUnigraph) <= 1) return(NULL)

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
	pathCheckLengths(from, to, arrows, values, free, labels, lbound, ubound, strictUnigraph, loop)

  numBounds <- max(length(lbound), length(ubound))
  for (bx in 1:numBounds) {
    l1 <- lbound[ 1 + (bx-1) %% length(lbound) ]
    u1 <- ubound[ 1 + (bx-1) %% length(ubound) ]
    if (!is.na(l1) && !is.na(u1) && l1 >= u1) {
      stop(paste("Lower bound", l1, "is greater than or equal to upper bound", u1))
    }
  }

	# create a new MxPath in the MxModel
	return(new("MxPath", unalteredFrom, unalteredTo, arrows, values, free, labels, lbound, ubound, connect, joinKey, step, strictUnigraph))
}

pathCheckToAndFrom <- function(from, to){
	# check for a missing to or from
	if (any(is.na(from)) || any(is.na(to))) {
		stop("The \'from\' field or the \'to\' field contains an NA", call. = FALSE)
	}
}

pathCheckLengths <- function(from, to, arrows, values,
        free, labels, lbound, ubound, strictUnigraph, loop) {
    numPaths <- max(length(from), length(to))
    pathCheckSingleLength(numPaths, length(arrows), "arrows", from, to, loop)
    pathCheckSingleLength(numPaths, length(values), "values", from, to, loop)
    pathCheckSingleLength(numPaths, length(free), "free/fixed designations", from, to, loop)
    pathCheckSingleLength(numPaths, length(labels), "labels", from, to, loop)
    pathCheckSingleLength(numPaths, length(lbound), "lbounds", from, to, loop)
    pathCheckSingleLength(numPaths, length(ubound), "ubounds", from, to, loop)
    pathCheckSingleLength(numPaths, length(strictUnigraph), "strictUnigraph", from, to, loop)
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
			deparse(width.cutoff = 400L, imxLocateFunction("mxPath")),
			omxQuotes(valname), "argument had class", class(value), "and length", length(value)
			), call. = FALSE)
	}
}

mxPath <- function(from, to = NA,
	connect = c("single", "all.pairs", "unique.pairs",
	            "all.bivariate", "unique.bivariate"), arrows = 1,
	free = TRUE, values = NA, labels = NA, lbound = NA, ubound = NA, ...,
  joinKey=as.character(NA), step=c(),
	strictUnigraph=TRUE) {
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
		msg <- paste("mxPath does not accept",
				omxQuotes(names(garbageArguments)),
				"as argument(s) and",
				"does not accept values",
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
  if (length(step)) {
    step <- as.integer(step)
    if (identical(arrows, 0) || (length(arrows) == length(step) &&
                                   all(arrows[!is.na(step)] == 0))) {
      # OK
    } else {
      stop(paste("step can only be provided for arrows=0 paths"))
    }
  }
  if (length(step)==0 && any(arrows==0)) {
    step <- rep(NA_integer_, length(arrows))
    step[arrows==0] <- 1L  # default step
  }
	if (isAllNa(to)) { to <- as.character(to) }
	if (isAllNa(from)) { from <- as.character(from) }
	if (isAllNa(values)) { values <- as.numeric(values) }
	if (isAllNa(labels)) { labels <- as.character(labels) }
	if (isAllNa(lbound)) { lbound <- as.numeric(lbound) }
	if (isAllNa(ubound)) { ubound <- as.numeric(ubound) }
	if (isAllNa(connect)) { connect <- as.character(connect) }
	pathCheckVector(from, 'from', is.character, 'character')
	pathCheckVector(to, 'to', is.character, 'character')
	pathCheckVector(arrows, 'arrows', is.numeric, 'numeric')
	pathCheckVector(free, 'free', is.logical, 'logical')
	pathCheckVector(labels, 'labels', is.character, 'character')
	if(!nchar0(labels)){
		stop("Found illegal label '', i.e. the empty label. Give us a real label, love.")
	}
	pathCheckVector(values, 'values', is.numeric, 'numeric')
	pathCheckVector(lbound, 'lbound', is.numeric, 'numeric')
	pathCheckVector(ubound, 'ubound', is.numeric, 'numeric')
	if (any(arrows==0)) pathCheckVector(step, 'step', is.numeric, 'numeric')
	pathCheckVector(strictUnigraph, 'strictUnigraph', is.logical, 'logical')
	generatePath(from, to, connect, arrows,
		values, free, labels,
		lbound, ubound, joinKey, step, strictUnigraph)
}

nchar0 <- function(x){
	all(nchar(x) > 0 | is.na(x))
}

prepPath <- function(path) {
	path@values[ is.na(path@values) ] <- 0

	if (single.na(path@to)) {
		# convert model.var -> var
		path@to <- sapply(path@from, function(x) {
			pieces <- strsplit(x, imxSeparatorChar, fixed = TRUE)[[1]]
			ifelse(length(pieces) == 2, pieces[2], pieces[1])
		}, USE.NAMES = FALSE)
	}

	expanded <- expandPathConnect(path@from, path@to, path@connect)
	path@from <- expanded$from
	path@to   <- expanded$to
	path
}

displayPath <- function(object) {
	cat(paste0("mxPath", '\n'))

	path <- prepPath(object)
	allfrom <- path@from
	allto <- path@to
	allarrows <- path@arrows
	allfree <- path@free
	allvalues <- path@values
	alllabels <- path@labels
	alllbound <- path@lbound
	allubound <- path@ubound
	allstep <- path@step
	#allstrictUnigraph <- path@strictUnigraph
	maxlength <- max(length(allfrom), length(allto))

	for(i in 0:(maxlength - 1)) {
		from <- allfrom[[i %% length(allfrom) + 1]]
		to <- allto[[i %% length(allto) + 1]]
		arrows <- allarrows[[i %% length(allarrows) + 1]]
		nextvalue <- allvalues[[i %% length(allvalues) + 1]]
		nextfree <- allfree[[i %% length(allfree) + 1]]
		nextlabel <- alllabels[[i %% length(alllabels) + 1]]
		nextubound <- allubound[[i %% length(allubound) + 1]]
		nextlbound <- alllbound[[i %% length(alllbound) + 1]]
    nextjoinkey <- path@joinKey[[i %% length(path@joinKey) + 1]]
    nextstep <- ifelse(length(allstep)==0, NA, allstep[[i %% length(allstep) + 1]])
    #nextstrictUnigraph <- allstrictUnigraph[[i %% length(alllabels) + 1]]

		cat(from)
		cat(paste0(' ', ifelse(arrows==1, "->", "<->"), ' '))
		cat(to)
		cat(paste0(" [value=",nextvalue))
		cat(paste0(", free=",nextfree))
		if (!is.na(nextlabel)) {
			cat(paste0(", label='", nextlabel,"'"))
		}
		if (!is.na(nextlbound)) {
			cat(paste0(", lbound=", nextlbound))
		}
		if (!is.na(nextubound)) {
			cat(paste0(", ubound=", nextubound))
		}
		if (!is.na(nextjoinkey)) {
			cat(paste0(", joinKey=", nextjoinkey))
		}
		if (!is.na(nextstep)) {
			cat(paste0(", step=", nextstep))
		}
		cat("]")
		cat('\n')
	}
}

setMethod("print", "MxPath", function(x,...) { displayPath(x) })
setMethod("show", "MxPath", function(object) { displayPath(object) })

#Maybe export this as an imx* function?:
AllRAMOrLISREL <- function(model,submodels=TRUE){
	out <- TRUE
	if(length(model$expectation)>0){
		# TODO: Add LISREL expectation once the backend knows how to calculate analytic matrix derivs for it:
		# out <- out && ( is(model$expectation,"MxExpectationRAM") || is(model$expectation,"MxExpectationLISREL") )
		out <- out && is(model$expectation,"MxExpectationRAM")
	}
	else{
		return(FALSE) #<--It has neither RAM (nor--TODO--LISREL) expectation, because it has no expectation.
	}
	if(submodels && length(model$submodels) > 0){
		out <- out && all(sapply(model@submodels, AllRAMOrLISREL))
	}
	return(out)
}

#Maybe export this as an imx* function?:
AnyRAMOrLISREL <- function(model,submodels=TRUE){
	out <- FALSE
	if(length(model$expectation)>0){
		# TODO: Add LISREL expectation once the backend knows how to calculate analytic matrix derivs for it:
		# out <- out || is(model$expectation,"MxExpectationRAM") || is(model$expectation,"MxExpectationLISREL")
		out <- out || is(model$expectation,"MxExpectationRAM")
	}
	else{
		if(!submodels){
			return(FALSE) #<--It has neither RAM (nor--TODO--LISREL) expectation, because it has no expectation.
		}
	}
	if(submodels && length(model$submodels) > 0){
		out <- out || all(sapply(model@submodels, AnyRAMOrLISREL))
	}
	return(out)
}

##' imxHasAlgebraOnPath
##'
##' This is an internal function exported for those people who know
##' what they are doing.  This function checks if a model (or any of its
##' submodels) either (1) has labels on MxPaths that reference one or more
##' MxAlgebras, or (2) defines any of the RAM matrices as MxAlgebras.
##'
##' @param model an MxModel object
##' @param submodels logical; recursion over child models?
##' @param strict logical; raise error if `model` contains no paths?
imxHasAlgebraOnPath <- function(model, submodels=TRUE, strict=FALSE){
	#out <- FALSE
	if(!AnyRAMOrLISREL(model,submodels=TRUE)){
		if(strict){
			#Throw error when strict, because model has no paths:
			stop(paste0(omxQuotes(model$name)," or one of its submodels does not use MxPaths")) 
		}
		else{
			return(FALSE) #<--Can't have algebras on paths if there are no paths!
		}
	}
	if(length(model@algebras)){
		allAlgNames <- names(model@algebras)
		#We need to check to see if the RAM matrices themselves are algebras:
		if( (model$expectation$A %in% allAlgNames) || (model$expectation$S %in% allAlgNames) || 
				(model$expectation$F %in% allAlgNames) ){
			#^^^Admittedly, it is difficult to imagine a scenario in which the user has specified the 'F' matrix as an *algebra*,
			#but we'd better check just in case.
			return(TRUE)
		}
		if( !is.na(model$expectation$M) && (model$expectation$M %in% allAlgNames) ){
			return(TRUE)
		}
		allPathLabels <- names(omxGetParameters(model,free=F)) #<--Paths that have labels referencing an algebra must be fixed.
		for(i in 1:length(allAlgNames)){
			#out <- out || length(grep(paste0(allAlgNames[i],"\\["),allPathLabels))
			#We are looking for the complete algebra name, followed by an opening square bracket:
			if(length(grep(paste0(allAlgNames[i],"\\["),allPathLabels))){ 
				return(TRUE)
			}
			# if(out){
			# 	return(out)
			# }
		}
	}
	if(submodels && length(model$submodels) > 0){
		return(any(sapply(model$submodels,imxHasAlgebraOnPath)))
	}
	return(FALSE) #<--Function couldn't find any reason to return TRUE.
}
