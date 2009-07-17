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
		model@options <- list()
		model@unsetoptions <- list()
		return(model)
	}
	model@options[[key]] <- value
	if (key %in% model@unsetoptions && !is.null(value)) {
		model@unsetoptions <- setdiff(model@unsetoptions, key)
	}
	if (is.null(value) && 
		key %in% names(getOption('mxOptimizerOptions'))) {
		model@unsetoptions <- union(model@unsetoptions, key)
	}
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
	"Linesearch tolerance" = "0.3",
	"Derivative Level" = "0",
	"Hessian" = "Yes"
)

# Convert the keys and values into strings
processOptionsList <- function(input) {
	keys <- sapply(names(input), as.character)
	values <- sapply(input, as.character)
	options <- list()
	options[keys] <- values
	return(options)
}

# Begin with the default NPSOL options,
# and then add the model-specific options.
generateOptionsList <- function(input, unset) {
	modified <- processOptionsList(input)
	options <- getOption('mxOptimizerOptions')
	options[as.character(unset)] <- NULL
	options[names(modified)] <- modified
	return(options)
}

