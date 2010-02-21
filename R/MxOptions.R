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

mxOption <- function(model, key, value, reset = FALSE) {
	if (reset) {
		model@options <- getOption('mxOptimizerOptions')
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
	"Calculate Hessian" = "No",
	"Standard Errors" = "No",
	"Chi-square Confidence Intervals" = "No"	
)

# Convert the keys and values into strings
generateOptionsList <- function(input) {
	options <- getOption('mxOptimizerOptions')
	if (length(input) > 0) {
		keys <- sapply(names(input), as.character)
		values <- sapply(input, as.character)
		options[keys] <- values
	}
	return(options)
}
