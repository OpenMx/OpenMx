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
	if (reset) {
		model@options <- list()
		return(model)
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
	"Function Precision" = "1e-14",
	"Infinite Bound Size" = "1.0e+15",
	"Feasibility tolerance" = "1.0e-05",
	"Major iterations" = "1000",
	"Verify level" = "3",
	"Line search tolerance" = "0.3",
	"Derivative Level" = "0",
	"Hessian" = "Yes",
	"Calculate Hessian" = "Yes",
	"Standard Errors" = "Yes",
	"CI Max Iterations" = "5"
)

checkpointOptions <- list(
	"Checkpoint Directory" = ".", 
	"Checkpoint Prefix" = "",
	"Checkpoint Units" = "minutes", 
	"Checkpoint Count Minutes" = 10,
	"Checkpoint Count Iterations" = 100,
	"Socket Server" = "", 
	"Socket Port" = 8080,
	"Socket Units" = "minutes", 
	"Socket Count Minutes" = 0.08,
	"Socket Count Iterations" = 1
)

otherOptions <- list(
	"No Sort Data" = character(),
	"RAM Optimization" = TRUE,
	"RAM Max Depth" = NA
)

generateOptionsList <- function(input, constraints) {
	if (is.null(input[["Standard Errors"]]) && length(constraints) > 0) {
		input[["Standard Errors"]] <- "No"
	}
	if (is.null(input[["Calculate Hessian"]]) && length(constraints) > 0) {
		input[["Calculate Hessian"]] <- "No"
	}
	options <- combineDefaultOptions(input)
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
