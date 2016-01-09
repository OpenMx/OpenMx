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

# returns a list object
# The names of the list are the full identifiers of named entities in the model.
#
# Given the ith element of the list, define the ith name as the source entity
# of the ith element.
#
# Each element of a list is a character vector.
# The character vector stores the named entities that are be affected by
# modifying the source entity.
transitiveClosure <- function(flatModel, dependencies) {
	dependencies <- dependencies@edges
	cache <- list()
    if (length(flatModel@matrices) == 0) {
		return(cache)
	}
	matrices <- names(flatModel@matrices)
	for (i in 1:length(matrices)) {
		target <- matrices[[i]]
		if (!(target %in% names(cache))) {
			cache <- transitiveClosureEntity(flatModel, dependencies, target, cache)
		}
	}
	return(cache)
}

extractElement <- function(name, object) { object[[name, exact=TRUE]] }

transitiveClosureEntity <- function(flatModel, dependencies, target, cache) {
	sinks <- dependencies[[target]]
	if (is.null(sinks)) {
		cache[[target]] <- character()
		return(cache)
	}
	isMissing <- !(sinks %in% names(cache))
	missing <- sinks[isMissing]
	if (length(missing) > 0) {
		for (i in 1:length(missing)) {
			entity <- missing[[i]]
			cache <- transitiveClosureEntity(flatModel, dependencies, entity, cache)
		}
	}
	entities <- lapply(sinks, extractElement, cache)
	combined <- c(entities, sinks)
	result <- Reduce(union, combined, character())
	cache[[target]] <- result
	return(cache)
}

doLocateIndex <- function(name, model, referant) {
	return(imxLocateIndex(model, name, referant))
}

generateDeps <- function(sinks, source) {
	retval <- rep.int(source, length(sinks))
	names(retval) <- sinks
	return(retval)
}

getPair <- function(val, x, y) {
	retval <- unlist(c(x[val], y[val]))
	names(retval) <- NULL
	return(retval)
}

combinePairs <- function(x, y) {
	shared <- intersect(names(x), names(y))
	diffx <- setdiff(names(x), names(y))
	diffy <- setdiff(names(y), names(x))
	common <- lapply(shared, getPair, x, y)
	names(common) <- shared
	xunique <- x[diffx]
	yunique <- y[diffy]
	return(c(common, xunique, yunique))
}

