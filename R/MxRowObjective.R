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

mxRowObjective <- function(rowAlgebra, reduceAlgebra, dimnames, rowResults = "rowResults", 
	filteredDataRow = "filteredDataRow", existenceVector = "existenceVector") {
	if (missing(rowAlgebra) || typeof(rowAlgebra) != "character") {
		stop("the 'rowAlgebra' argument is not a string (the name of the row-by-row algebra)")
	}
	if (missing(reduceAlgebra) || typeof(reduceAlgebra) != "character") {
		stop("the 'reduceAlgebra' argument is not a string (the name of the reduction algebra)")
	}
	if (missing(dimnames) || typeof(dimnames) != "character") {
		stop("the 'dimnames' argument is not a string (the column names from the data set)")
	}
	if (any(is.na(dimnames))) {
		stop("NA values are not allowed for 'dimnames' vector")
	}
	rowResults <- checkStringArgument(rowResults, "rowResults")
	filteredDataRow <- checkStringArgument(filteredDataRow, "filteredDataRow")
	existenceVector <- checkStringArgument(existenceVector, "existenceVector")
	expectation <- NULL	
	fitfunction <- mxFitFunctionRow(rowAlgebra, reduceAlgebra, dimnames, 
		rowResults, filteredDataRow, existenceVector)
	msg <- paste("Objective functions have been deprecated.",
		"Please use mxFitFunctionRow() instead.")
	warning(msg)
	return(list(expectation=expectation, fitfunction=fitfunction))
}


