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

mxFIMLObjective <- function(covariance, means, dimnames = NA, 
	thresholds = NA, vector = FALSE, threshnames = dimnames) {
	if (missing(covariance) || typeof(covariance) != "character") {
		stop("'covariance' argument is not a string (the name of the expected covariance matrix)")
	}
	if (missing(means) || typeof(means) != "character") {
		stop("'means' argument is not a string (the name of the expected means vector)")
	}
	if (single.na(thresholds)) thresholds <- as.character(NA)
	if (single.na(dimnames)) dimnames <- as.character(NA)
	if (!is.vector(dimnames) || typeof(dimnames) != 'character') {
		stop("'dimnames' argument is not a character vector")
	}
	if (length(thresholds) != 1) {
		stop("'thresholds' argument must be a single matrix or algebra name")
	}
	if (length(dimnames) == 0) {
		stop("'dimnames' argument cannot be an empty vector")
	}
	if (length(dimnames) > 1 && any(is.na(dimnames))) {
		stop("NA values are not allowed for 'dimnames' vector")
	}
	if (length(vector) > 1 || typeof(vector) != "logical") {
		stop("'vector' argument is not a logical value")
	}
	threshnames <- checkThreshnames(threshnames)
	expectation <- mxExpectationNormal(covariance, means, dimnames, 
		thresholds, threshnames)
	fitfunction <- mxFitFunctionML(vector)
	msg <- paste("Objective functions like mxFIMLObjective() have been deprecated in favor of expectation and fit functions.\n",
		"Please use mxExpectationNormal(covariance= , means = , ...) instead, and add a call to mxFitFunctionML(). See examples at help(mxExpectationNormal)")
	warning(msg)
	return(list(expectation=expectation, fitfunction=fitfunction))
}

