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

# returns a list of paths
generatePath <- function(from, to, 
		all, arrows, values, free,
		labels, lbound, ubound) {
		if (single.na(to)) {
			to <- from
		}
		from <- as.list(from)
		to <- as.list(to)
		if (all) {
			from <- rep(from, each=length(to))	
		}
		result <- suppressWarnings(mapply(generateSinglePath, from, to, 
		arrows, values, free,
		labels, lbound, ubound, SIMPLIFY = FALSE))
		return(result)
}

generateSinglePath <- function(from, to, 
		arrows, values, free,
		labels, lbound, ubound) {
	result <- list()
	result[['from']] <- from
	result[['to']] <- to
	result[['arrows']] <- arrows[[1]]
	result[['values']] <- values[[1]]	
	result[['free']] <- free[[1]]
	result[['labels']] <- labels[[1]]	
	result[['lbound']] <- lbound[[1]]
	result[['ubound']] <- ubound[[1]]
	return(result)
}

omxIsPath <- function(value) {
	return(is.list(value) && 
		!is.null(value[['from']]) &&
		!is.null(value[['to']]))
}


mxPath <- function(from, to = NA, all = FALSE, arrows = 1, free = TRUE,
	values = NA, labels = NA, lbound = NA, ubound = NA) {
	if (length(values) == 1 && is.na(values)) values <- NULL
	if (length(labels) == 1 && is.na(labels)) labels <- NULL
	if (length(lbound) == 1 && is.na(lbound)) lbound <- NULL
	if (length(ubound) == 1 && is.na(ubound)) ubound <- NULL
	generatePath(from, to, all, arrows, 
		values, free, labels, 
		lbound, ubound)
}
