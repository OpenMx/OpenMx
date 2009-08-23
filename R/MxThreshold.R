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

convertThresholds <- function(observed, thresholds, objname) {
	if (dimnames(thresholds) == NULL || dimnames(thresholds)[[2]] == NULL) {
		stop(paste("The thresholds matrix for objective function", 
			omxQuotes(objname), "does not contain column names"), call. = FALSE)
	}
	if (dimnames(observed) == NULL || dimnames(observed)[[2]] == NULL) {
		stop(paste("The observed data matrix for objective function", 
			omxQuotes(objname), "does not contain column names"), call. = FALSE)
	}
	threshNames <- dimnames(observed)[[2]]
	obsNames <- dimnames(observed)[[2]]
	missingNames <- setdiff(threshNames, obsNames)
	if (length(missingNames) > 0) {
		stop(paste("The following column names in the threshold",
			"matrix do not exist in the observed data matrix",
			"for objective function", omxQuotes(objname), 
			":", omxQuotes(missingNames)), call. = FALSE)
	}
	for(i in 1:length(threshNames)) {
		tColumn <- as.ordered(observed[,threshNames[[i]]])
		observed[,threshNames[[i]]] <- tColumn
		tCount <- sum(!is.na(tColumn))
		if (tCount != (length(levels(tColumn)) - 1)) {
			stop(paste("The number of thresholds in column",
				omxQuotes(threshNames[[i]]),
				"is not one less than the number of levels",
				"in objective function", 
				omxQuotes(objname)), call. = FALSE)
		}
		tValues <- tColumn[1:tCount]
		if (any(is.na(tValues))) {
			stop(paste("The thresholds in column",
				omxQuotes(threshNames[[i]]),
				"contain NA values in between non-NA values",
				"in objective function",
				omxQuotes(objname)), call. = FALSE)			
		}
		if (tCount < length(tColumn) && 
				any(!is.na(tColumn[tCount + 1:length(tColumn)]))) {
			stop(paste("The thresholds in column",
				omxQuotes(threshNames[[i]]),
				"contain NA values in between non-NA values",
				"in objective function",
				omxQuotes(objname)), call. = FALSE)
		}
		sortValues <- sort(tValues)
		if (!all(sortValues == tValues)) {
			stop(paste("The thresholds in column",
				omxQuotes(threshNames[[i]]),
				"are not in sorted order",
				"in objective function",
				omxQuotes(objname)), call. = FALSE)	
		}
	}
	return(observed)
}