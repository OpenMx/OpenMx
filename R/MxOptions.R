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

mxOption <- function(model, key, value, reset = FALSE) {
	if (length(model) != 1 || !is(model, "MxModel")) {
		stop("argument 'model' must be an MxModel object")
	}
	if (length(reset) != 1 || !is.logical(reset)) {
		stop("argument 'reset' must be TRUE or FALSE")
	}
	if (reset) {
		model@options <- list()
		return(model)
	}
	if (length(key) != 1 || !is.character(key)) {
		stop("argument 'key' must be a character string")
	}
	if (length(value) > 1) {
		stop("argument 'value' must be either NULL or of length 1")
	}
	optionsNames <- names(getOption('mxOptions'))
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
			"the capitalization."))
	}
	model@options[[key]] <- value
	return(model)
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
	"Major iterations" = "1000",
	"Verify level" = "3",
	"Line search tolerance" = "0.3",
	"Derivative level" = "0",
	"Hessian" = "Yes",
	"Calculate Hessian" = "Yes",
	"Standard Errors" = "Yes",
	"CI Max Iterations" = "5",
	"Number of Cores" = 1
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
	"Error Checking" = "Yes",
	"No Sort Data" = character(),
	"RAM Inverse Optimization" = "Yes",
	"RAM Max Depth" = NA
)

generateOptionsList <- function(input, constraints, useOptimizer) {
	if (is.null(input[["Standard Errors"]]) && length(constraints) > 0) {
		input[["Standard Errors"]] <- "No"
	}
	if (is.null(input[["Calculate Hessian"]]) && length(constraints) > 0) {
		input[["Calculate Hessian"]] <- "No"
	}
	options <- combineDefaultOptions(input)
	if (useOptimizer) {
		options[["useOptimizer"]] <- "Yes"
	} else {
		options[["useOptimizer"]] <- "No"
	}
	return(options)
}

# Convert the keys and values into strings
combineDefaultOptions <- function(input) {
	options <- getOption('mxOptions')
	options <- options[names(options) %in% names(npsolOptions)]
	input <- input[names(input) %in% names(npsolOptions)]
	if (length(input) > 0) {
		keys <- sapply(names(input), as.character)
		values <- sapply(input, as.character)
		options[keys] <- values
	}
	return(options)
}
