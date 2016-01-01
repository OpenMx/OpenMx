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

mxRObjective <- function(objfun, ...) {
	if (!is.function(objfun)) {
		stop("First argument 'objfun' must be of type function")
	}
	if (length(formals(objfun)) != 2) {
		stop("The objective function must take exactly two arguments: a model and a persistant state")
	}
	expectation <- NULL
	fitfunction <- mxFitFunctionR(objfun, ...)
	msg <- paste("Objective functions have been deprecated.",
		"Please use mxFitFunctionR() instead.")
	warning(msg)
	return(list(expectation=expectation, fitfunction=fitfunction))
}


