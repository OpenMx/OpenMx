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

mxRAMObjective <- function(A, S, F, M = NA, dimnames = NA, thresholds = NA, vector = FALSE,
									threshnames = dimnames) {
	if (missing(A) || typeof(A) != "character") {
		msg <- paste("argument 'A' is not a string",
			"(the name of the 'A' matrix)")
		stop(msg)
	}	
	if (missing(S) || typeof(S) != "character") {
		msg <- paste("argument 'S' is not a string",
			"(the name of the 'S' matrix)")
		stop(msg)
	}
	if (missing(F) || typeof(F) != "character") {
		msg <- paste("argument 'F' is not a string",
			"(the name of the 'F' matrix)")
		stop(msg)
	}
	if (!(single.na(M) || typeof(M) == "character")) {
		msg <- paste("argument M is not a string",
			"(the name of the 'M' matrix)")
		stop(msg)
	}
	if (is.na(M)) M <- as.integer(NA)
	if (single.na(thresholds)) thresholds <- as.character(NA)
	if (single.na(dimnames)) dimnames <- as.character(NA)
	if (!is.vector(dimnames) || typeof(dimnames) != 'character') {
		stop("Dimnames argument is not a character vector")
	}
	if (length(thresholds) != 1) {
		stop("Thresholds argument must be a single matrix or algebra name")
	}
	if (length(dimnames) == 0) {
		stop("Dimnames argument cannot be an empty vector")
	}
	if (length(dimnames) > 1 && any(is.na(dimnames))) {
		stop("NA values are not allowed for dimnames vector")
	}
	if (length(vector) > 1 || typeof(vector) != "logical") {
		stop("Vector argument is not a logical value")
	}
	threshnames <- checkThreshnames(threshnames)
	expectation <- mxExpectationRAM(A, S, F, M, dimnames, thresholds, threshnames)
	fitfunction <- mxFitFunctionML(vector)
	msg <- paste("Objective functions like mxRAMObjective() have been deprecated in favor of expectation and fit functions.\n",
		"Please use mxExpectationRAM(A, S, F, M, ...) instead, and add a call to mxFitFunctionML().\n",
		"See examples at help(mxExpectationRAM)")
	warning(msg)
	return(list(expectation=expectation, fitfunction=fitfunction))
}


