#
#   Copyright 2007-2009 The OpenMx Project
#
#   Licensed under the Apache License, Version 2.0 (the "License");
#   you may not use this file except in compliance with the License.
#   You may obtain a copy of the License at
# 
#        http://www.apache.org/licenses/LICENSE-2.0
# 
#    Unless required by applicable law or agreed to in writing, software
#    distributed under the License is distributed on an "AS IS" BASIS,
#    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#   See the License for the specific language governing permissions and
#   limitations under the License.

# returns a list of paths
omxPath <- function(from, to = NULL, all = FALSE, free = NULL, 
	arrows = NULL, startVal = NULL, 
	endVal = NULL, algebra = NULL,
	name = NULL, label = NULL,
	boundMax = NULL, boundMin = NULL,
	ciUpper = NULL, ciLower = NULL) {
		if (is.null(to)) {
			to <- from
		}
		from <- as.list(from)
		to <- as.list(to)
		if (all) {
			from <- rep(from, each=length(to))	
		}
		result <- suppressWarnings(mapply(omxSinglePath, from, to,
			free, arrows, startVal, endVal,
				algebra, name, label, boundMax,
				boundMin, ciUpper, ciLower, SIMPLIFY=FALSE))
		return(result)
}

omxSinglePath <- function(from, to, free = NULL, 
	arrows = NULL, startVal = NULL, 
	endVal = NULL, algebra = NULL,
	name = NULL, label = NULL,
	boundMax = NULL, boundMin = NULL,
	ciUpper = NULL, ciLower = NULL) {
	result <- list()
	result[['from']] <- from
	result[['to']] <- to
	result[['free']] <- free[[1]]
	result[['arrows']] <- arrows[[1]]
	result[['startVal']] <- startVal[[1]]
	result[['endVal']] <- endVal[[1]]
	result[['algebra']] <- algebra[[1]]
	result[['name']] <- name[[1]]
	result[['label']] <- label[[1]]	
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


mxPath <- function(from, to = NULL, all = FALSE, free = NULL, 
	arrows = NULL, startVal = NULL, 
	endVal = NULL, algebra = NULL,
	name = NULL, label = NULL,
	boundMax = NULL, boundMin = NULL,
	ciUpper = NULL, ciLower = NULL) {

	omxPath(from, to, all, free, 
		arrows, startVal, endVal, 
		algebra, name, label, boundMax,
		boundMin, ciUpper, ciLower)
}
