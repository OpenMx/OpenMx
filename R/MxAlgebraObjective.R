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

mxAlgebraObjective <- function(algebra, numObs = NA, numStats = NA) {
	if (missing(algebra) || typeof(algebra) != "character") {
		stop("Algebra argument is not a string (the name of the algebra)")
	}
	if (single.na(numObs)) {
		numObs <- as.numeric(NA)
	}
	if (single.na(numStats)) {
		numStats <- as.numeric(NA)
	}
	expectation <- NULL
	fitfunction <- mxFitFunctionAlgebra(algebra, numObs, numStats)
	msg <- paste("Objective functions like mxAlgebraObjective() have been deprecated in favor of expectation and fit functions.\n",
		"Please use mxFitFunctionAlgebra(algebra = ...). See examples at help(mxFitFunctionAlgebra)")
	warning(msg)
	return(list(expectation = expectation, fitfunction = fitfunction))
}
