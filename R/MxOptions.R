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

mxOption <- function(model=NULL, key=NULL, value, reset = FALSE) {
	if (!reset && (length(key) != 1 || !is.character(key))) {
		stop("argument 'key' must be a character string")
	}
	if (!missing(model) && !is.null(model) && !is(model, "MxModel")) {
		stop(paste("The first argument to mxOption must",
			"be an MxModel, not", omxQuotes(class(model))))
	}
	if (is.null(model) && reset) {
		return(invisible(mxSetDefaultOptions()))
	}
	if (missing(value)) {
		if (length(model) && !is.null(model@options[[key]])) {
			return(model@options[[key]])
		}
		return(processDefaultOptionList(key, value))
	}
	if (length(value) > 1 && key!="No Sort Data" && key != "Status OK") {
		msg <- paste("argument 'value' must be either NULL or of length 1.",
			"You gave me an object of length", length(value))
		stop(msg)
	}
	if (length(reset) != 1 || !is.logical(reset)) {
		stop("argument 'reset' must be TRUE or FALSE")
	}
	if (key == "Major iterations" && typeof(value) == "closure") {
		args <- formals(value)
		if (length(args) != 2) {
			msg <- paste("The function provided to the option 'Major iterations'",
				"must have exactly 2 arguments but you have provided",
				"a function with", length(args), "arguments.")
			stop(msg)
		}
		if (!single.na(match("...", names(args)))) {
			msg <- paste("You have provided a function to the option 'Major iterations'",
				"that uses the '...' argument.")
			stop(msg)
		}
	}
    if (is.null(model)) {
	    return(processDefaultOptionList(key, value))
    }
	if (length(model) > 1 || !is(model, "MxModel")) {
		stop("argument 'model' must be an MxModel object")
	}
	if (reset) {
		model@options <- list()
		return(model)
	}
	optionsNames <- names(getOption('mxOptions'))
	match <- grep(paste("^", key, "$", sep = ""), optionsNames,
		ignore.case=TRUE)
	if(length(match) == 0) {
		stop(paste("argument 'key' is the character string",
			omxQuotes(key), "and cannot be found in",
			"getOption('mxOptions')"))
	}
	key <- optionsNames[[match]] # repair capitalization
	if (key == "Default optimizer" || key == "Gradient algorithm" || key == "Gradient iterations") {
		stop(paste(omxQuotes(key), " is a global option and cannot be set on models.\n",
		"To change ", omxQuotes(key) ," globally, use, e.g.:\n",
		"mxOption(NULL, '", key, "', '", value,"')", sep = ""))
        # to use NLOPT, use: mxOption(NULL, 'Default optimizer', 'NLOPT')
	}
	if (key == "Status OK") value <- as.statusCode(value)
	model@options[[key]] <- value
	return(model)
}

processDefaultOptionList <- function(key, value) {
	defaultOptions <- getOption('mxOptions')
	optionsNames <- names(defaultOptions)
	match <- grep(paste("^", key, "$", sep = ""), optionsNames,
		ignore.case=TRUE)
	if(length(match) == 0) {
		stop(paste("argument 'key' has a value",
			omxQuotes(key), "that cannot be found in",
			"getOption('mxOptions')"))
	}
	key <- optionsNames[[match]] # repair capitalization
	if (missing(value)) return(defaultOptions[[key]])
	defaultOptions[[key]] <- value
	options('mxOptions' = defaultOptions)
	return(invisible(defaultOptions))
}

##' imxDetermineDefaultOptimizer
##'
##' This is an internal function exported for those people who know
##' what they are doing.
##'
##' @details Returns a character, the default optimizer
imxDetermineDefaultOptimizer <- function() {
	engine <- Sys.getenv("IMX_OPT_ENGINE")
	if (!nchar(engine)) {
		if (imxHasNPSOL()) {
			engine <- "SLSQP"
		} else {
			engine <- "SLSQP"
		}
	}
	engine
}

# Names and values must all be strings
npsolOptions <- list(
	"Nolist" = "",
	"Print level" = "0",
	"Minor print level" = "0",
	"Print file" = "0",
	"Summary file" = "0",
	"Function precision" = "Auto",#"1e-14"
	"Optimality tolerance" = "6.3e-12",
	"Infinite bound size" = "1.0e+15",
	"Feasibility tolerance" = "5e-2",
	"Major iterations" = function(nParams, nConstraints) { max(1000, 3 * nParams + 10 * nConstraints) },
	"Verify level" = "-1",
	"Line search tolerance" = "0.3",
	"Derivative level" = "0",
    "Step limit" = "2.0",
	"Hessian" = "Yes",
# below are not npsol options
	"Calculate Hessian" = "Yes",
	"Standard Errors" = "Yes",
	"Analytic Gradients" = "Yes"
)

checkpointOptions <- list(
	"Checkpoint Directory" = ".", 
	"Checkpoint Prefix" = "",
	"Checkpoint Units" = "iterations",
	"Checkpoint Count" = 1,
        "Checkpoint Fullpath" = "",
	"Socket Server" = "", 
	"Socket Port" = 8080,
	"Socket Units" = "minutes", 
	"Socket Count" = c("minutes" = 0.08, "iterations" = 1)
)

otherOptions <- list(
    "Always Checkpoint" = "No",
	"Error Checking" = "Yes",
	"No Sort Data" = character(),
	"RAM Inverse Optimization" = "Yes",
	"RAM Max Depth" = NA,
	"UsePPML" = "No",
	"Allow Unlabeled" = FALSE,
    "loglikelihoodScale" = -2.0,
    "maxOrdinalPerBlock" = 20,
    "mvnMaxPointsA" = 0,
    "mvnMaxPointsB" = 0,
    "mvnMaxPointsC" = 0,
    "mvnMaxPointsD" = 3.606464,
    "mvnMaxPointsE" = -0.126859,
    "mvnAbsEps" = 0,
    "mvnRelEps" = .005,
    "maxStackDepth" = 25000L,   # R_PPSSIZE/2
    "Gradient algorithm" = NULL,
    "Gradient iterations" = "Auto",#1L,
    "Gradient step size" = "Auto",#1.0e-7,
    "Parallel diagnostics" = "No",
    "Debug protect stack" = "No",
    "Nudge zero starts" = "Yes",
    "Status OK"= as.statusCode(c("OK", "OK/green")),
    "Max minutes"=0
)

limitMajorIterations <- function(options, numParam, numConstraints) {
	mIters <- options[["Major iterations"]]
	if (typeof(mIters) == "closure") {
		mIters <- do.call(mIters, list(numParam, numConstraints))
	}
	options[["Major iterations"]] <- as.character(mIters)
	options
}

imxGetNumThreads <- function() {
	if (imxSfClient()) {
		return(1L)
	} else {
		thrlimit <- as.integer(Sys.getenv("OMP_NUM_THREADS"))
		if (!is.na(thrlimit)) {
			return(thrlimit)
		} else {
			detect <- omxDetectCores()
			if(is.na(detect)) detect <- 1L
					# Due to demand by CRAN maintainers, we default to 2 cores
					# when OMP_NUM_THREADS is not set. This seems like a bad
					# policy to the OpenMx team, but we have no choice.
			else detect <- 2L
			return(detect)
		}
	}
}

generateOptionsList <- function(model, constraints, useOptimizer) {
	input <- list()
	if (!is.null(model)) {
		input <- model@options
		if (is.null(input[["Standard Errors"]]) && constraints > 0 && imxHasWLS(model)){
			input[["Standard Errors"]] <- "No"
		}
		if( !is.null(input[["UsePPML"]]) 
		   && (input[["UsePPML"]] == "PartialSolved" || input[["UsePPML"]] == "Split") ) {
			input[["Calculate Hessian"]] <- "No"
			input[["Hessian"]] <- "No"
			input[["Standard Errors"]] <- "No"
		}
	}
	options <- combineDefaultOptions(input)
	if (useOptimizer) {
		options[["useOptimizer"]] <- "Yes"
		#PPML Analytical solution
		if (!is.null(model@options$UsePPML) && model@options$UsePPML == "Solved")
			options[["useOptimizer"]] <- "No"
	} else {
		options[["useOptimizer"]] <- "No"
	}
	if (identical(options[["Standard Errors"]], "Yes") &&
		identical(options[["Calculate Hessian"]], "No")) {
		msg <- paste('The "Standard Errors" option is enabled and',
		'the "Calculate Hessian" option is disabled. This may',
			     'result in poor accuracy standard errors.')
		warning(msg)
	}
	return(options)
}

# Convert the keys and values into strings
combineDefaultOptions <- function(input) {
	options <- getOption('mxOptions')
	temp <- input[names(input) %in% names(npsolOptions)]
	temp[["Major iterations"]] <- NULL
	if (length(temp) > 0) {
		keys <- sapply(names(temp), as.character)
		values <- sapply(temp, as.character)
		ynOptions <- options[keys]=='Yes' | options[keys]=='No'
		badYN <- values[ynOptions] != 'Yes' & values[ynOptions] != 'No'
		if (any(badYN)) {
			stop(paste("mxOption '", names(badYN),
				   "' must be either 'Yes' or 'No'\n", sep=''))
		}
		options[keys] <- values
	}
	if (!is.null(input[["Major iterations"]])) {
		options[["Major iterations"]] <- input[["Major iterations"]]
	}
  #Need to make sure that non-default values for options not already handled in this function don't get
  #overwritten by the defaults:
  namesHandled <- c( names(temp), "Major iterations" )
	if(sum( !(names(input) %in% namesHandled) )>0){
    options[names(input)[!(names(input) %in% namesHandled)]] <- 
      input[names(input)[!(names(input) %in% namesHandled)]]
  }
	return(options)
}


##' imxAutoOptionValue
##' 
##' Convert "Auto" placeholders in global mxOptions to actual default values.
##' 
##' This is an internal function exported for documentation purposes.
##' Its primary purpose is to convert the on-load value of "Auto"to
##' valid values for \link{mxOption}s \sQuote{Gradient step size},
##' \sQuote{Gradient iterations}, and 
##' \sQuote{Function precision}--respectively, 1.0e-7, 1L, and 1e-14.
##' 
##' @param optionName Character string naming the \link{mxOption} for which a numeric or integer value is wanted.
##' @param optionList List of options; defaults to list of global \link{mxOption}s.
##' imxAutoOptionValue
imxAutoOptionValue <- function(optionName, optionList=options()$mxOption){
	#First, check to see if the option already has a valid value (possibly in string form), and if so, return that:
	numcast <- try(suppressWarnings(as.numeric(optionList[[optionName]])),silent=TRUE)
	if(!length(numcast)){
		#NULL numcast is most likely to result from either (1) misspelled optionName, 
		#or (2) user providing non-default value for optionList that lacks an element named optionName.
		#Throwing an error seems the best behavior in this case.
		stop(paste("extracting element '",optionName,"' from argument 'optionList' resulted in NULL"),sep="")
	}
	#numcast will be try-error for e.g. NPSOL option "Major iterations" (on-load default is a function);
	if("try-error" %in% class(numcast)){return(optionList[[optionName]])}
	if(length(numcast) && !is.na(numcast)){
		if(optionName=="Gradient iterations"){numcast <- as.integer(numcast)}
		return(numcast)
	}
	#Otherwise, if the current value is a string and can be matched to "Auto",
	#convert to default numerical value for the three motivating cases: 
	else{
		if(length(grep(pattern="Auto",x=optionList[[optionName]],ignore.case=T))){
			if(optionName=="Gradient step size"){return(1.0e-7)}
			if(optionName=="Gradient iterations"){return(1L)}
			if(optionName=="Function precision"){return(1e-14)}
		}
		else{stop(paste("found unrecognized character string '",optionList[[optionName]],"' as value for mxOption '",optionName,"' in argument 'optionList'",sep=""))}
	}
}

