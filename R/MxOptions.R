#
#   Copyright 2007-2014 The OpenMx Project
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

mxOption <- function(model, key, value, reset = FALSE) {
	if (length(key) != 1 || !is.character(key)) {
		stop("argument 'key' must be a character string")
	}
	if (length(value) > 1 && key!="No Sort Data") {
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
    if (length(model) == 0 && is.null(model)) {
        return(processDefaultOptionList(key, value, reset))
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
	if (!identical(optionsNames[[match]], key)) {
		stop(paste("argument 'key' is the character string",
			omxQuotes(key), "but the option is named",
			omxQuotes(optionsNames[[match]]), ": please correct",
			"the capitalization and re-run mxOption()."))
	}
	model@options[[key]] <- value
	return(model)
}

processDefaultOptionList <- function(key, value, reset) {
	defaultOptions <- getOption('mxOptions')
	optionsNames <- names(defaultOptions)
	match <- grep(paste("^", key, "$", sep = ""), optionsNames,
		ignore.case=TRUE)
	if(length(match) == 0) {
		stop(paste("argument 'key' has a value",
			omxQuotes(key), "that cannot be found in",
			"getOption('mxOptions')"))
	}
	if (!identical(optionsNames[[match]], key)) {
		stop(paste("argument 'key' has a value",
			omxQuotes(key), "but the option is named",
			omxQuotes(optionsNames[[match]]), ": please correct",
			"the capitalization and re-run mxOption()."))
	}
	defaultOptions[[key]] <- value
	options('mxOptions' = defaultOptions)
	return(invisible(defaultOptions))
}

determineDefaultOptimizer <- function() {
	engine <- Sys.getenv("IMX_OPT_ENGINE")
	if (!nchar(engine)) {
		if (imxHasNPSOL()) {
			engine <- "NPSOL"
		} else {
			engine <- "CSOLNP"
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
	"Function precision" = "1e-14",
	"Optimality tolerance" = "6.3e-12",
	"Infinite bound size" = "1.0e+15",
	"Feasibility tolerance" = "1.0e-05",
	"Optimality tolerance" = as.character(1e-14 ^ 0.8),
	"Major iterations" = function(nParams, nConstraints) { max(1000, 3 * nParams + 10 * nConstraints) },
	"Verify level" = "-1",
	"Line search tolerance" = "0.3",
	"Derivative level" = "0",
	"Hessian" = "Yes",
# below are not npsol options
	"Calculate Hessian" = "Yes",
	"Standard Errors" = "Yes",
	"CI Max Iterations" = "5",
	"Analytic Gradients" = "Yes",
	"Number of Threads" = 0
)

checkpointOptions <- list(
	"Checkpoint Directory" = ".", 
	"Checkpoint Prefix" = "",
	"Checkpoint Units" = "minutes", 
	"Checkpoint Count" = c("minutes" = 10, "iterations" = 100),
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
    "maxptsa" = 1.1,
    "maxptsb" = 1.2,
    "maxptsc" = 1.3,
    "maxStackDepth" = 25000L   # R_PPSSIZE/2
)

generateOptionsList <- function(model, numParam, constraints, useOptimizer) {
	input <- model@options
	if (is.null(input[["Standard Errors"]]) && length(constraints) > 0) {
		input[["Standard Errors"]] <- "No"
	}
	if (is.null(input[["Calculate Hessian"]]) && length(constraints) > 0) {
		input[["Calculate Hessian"]] <- "No"
	}
	if( !is.null(input[["UsePPML"]]) 
		&& (input[["UsePPML"]] == "PartialSolved" || input[["UsePPML"]] == "Split") ) {
		input[["Calculate Hessian"]] <- "No"
		input[["Hessian"]] <- "No"
		input[["Standard Errors"]] <- "No"
	}
	options <- combineDefaultOptions(input)
	mIters <- options[["Major iterations"]]
	if (typeof(mIters) == "closure") {
		mIters <- do.call(mIters, list(numParam, length(constraints)))
	}
	options[["Major iterations"]] <- as.character(mIters)
	if (useOptimizer) {
		options[["useOptimizer"]] <- "Yes"
		#PPML Analytical solution
		if (!is.null(model@options$UsePPML) && model@options$UsePPML == "Solved")
			options[["useOptimizer"]] <- "No"
	} else {
		options[["useOptimizer"]] <- "No"
	}
	if (model@.forcesequential) {
		options[["Number of Threads"]] <- 1L 
	} else if (is.null(options[["Number of Threads"]]) || 
			options[["Number of Threads"]] == 0) {
		if (imxSfClient()) {
 			options[["Number of Threads"]] <- 1L 
		} else {
			thrlimit <- as.integer(Sys.getenv("OMP_NUM_THREADS"))
			if (!is.na(thrlimit)) {
				options[["Number of Threads"]] <- thrlimit
			} else {
				detect <- omxDetectCores()
				if(is.na(detect)) detect <- 1L
				options[["Number of Threads"]] <- detect 
			}
		}
	}
	if (identical(options[["Standard Errors"]], "Yes") &&
		identical(options[["Calculate Hessian"]], "No")) {
		msg <- paste('The "Standard Errors" option is enabled and',
		'the "Calculate Hessian" option is disabled. Generating',
		'standard errors requires the Hessian calculation. Please',
		'disable standard errors or enable the Hessian calculation.')
		stop(msg, call. = FALSE)
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
	return(options)
}
