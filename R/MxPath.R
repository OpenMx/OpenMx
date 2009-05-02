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
omxPath <- function(from, to, all, free, 
	arrows, start, name, description,
	boundMax, boundMin, ciUpper, ciLower) {
		if (single.na(to)) {
			to <- from
		}
		from <- as.list(from)
		to <- as.list(to)
		if (all) {
			from <- rep(from, each=length(to))	
		}
		result <- suppressWarnings(mapply(omxSinglePath, from, to,
					free, arrows, start, 
					name, description, boundMax,
					boundMin, ciUpper, ciLower, SIMPLIFY = FALSE))
		return(result)
}

omxSinglePath <- function(from, to, free, 
	arrows, start, name, description, 
	boundMax, boundMin,
	ciUpper, ciLower) {
	result <- list()
	result[['from']] <- from
	result[['to']] <- to
	result[['free']] <- free[[1]]
	result[['arrows']] <- arrows[[1]]
	result[['start']] <- start[[1]]
	result[['name']] <- name[[1]]
	result[['description']] <- description[[1]]	
	result[['boundMax']] <- boundMax[[1]]
	result[['boundMin']] <- boundMin[[1]]
	result[['ciUpper']] <- ciUpper[[1]]
	result[['ciLower']] <- ciLower[[1]]
	return(result)
}

omxIsPath <- function(value) {
	return(is.list(value) && 
		!is.null(value[['from']]) &&
		!is.null(value[['to']]))
}


mxPath <- function(from, to = NA, all = FALSE, free = TRUE, 
	arrows = 1, start = NA, name = NA, description = NA, 
	boundMax = NA, boundMin = NA,
	ciUpper = NA, ciLower = NA) {
	if (length(start) == 1 && is.na(start)) start <- NULL
	if (length(name) == 1 && is.na(name)) name <- NULL
	if (length(description) == 1 && is.na(description)) description <- NULL
	if (length(boundMax) == 1 && is.na(boundMax)) boundMax <- NULL
	if (length(boundMin) == 1 && is.na(boundMin)) boundMin <- NULL
	if (length(ciUpper) == 1 && is.na(ciUpper)) ciUpper <- NULL
	if (length(ciLower) == 1 && is.na(ciLower)) ciLower <- NULL
	omxPath(from, to, all, free, 
		arrows, start, name, 
		description, boundMax,
		boundMin, ciUpper, ciLower)
}
