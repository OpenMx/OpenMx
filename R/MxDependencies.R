#
#   Copyright 2007-2012 The OpenMx Project
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

doLocateIndex <- function(name, model, referant) {
	return(imxLocateIndex(model, name, referant))
}

entityDependencies <- function(entity, flatModel) {
	dependencies <- names(entity@dependencies)
	if (length(dependencies) > 0) {
		dependencies <- sapply(dependencies, doLocateIndex, flatModel, entity@name) 
	} else {
		dependencies <- integer(0)
	}
	entity@dependencies <- dependencies
	return(entity)
}

convertDependencies <- function(flatModel) {
	flatModel@algebras <- lapply(flatModel@algebras, entityDependencies, flatModel)
	flatModel@objectives <- lapply(flatModel@objectives, entityDependencies, flatModel)
	return(flatModel)
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

applyDeps <- function(object, pairs) {
	sources <- pairs[[object@name]]
	if (is.null(sources)) {
		return(object)
	} else {
		deplist <- rep.int(1L, length(sources))
		names(deplist) <- sources
		object@dependencies <- deplist
	}
	return(object)
}

addDependencies <- function(flatModel, dependencies) {
	edges <- dependencies@edges
	if (length(edges) == 0) {
		return(flatModel)
	}
	pairs <- mapply(generateDeps, edges, names(edges), USE.NAMES = FALSE)
	pairs <- Reduce(combinePairs, pairs, character())
	flatModel@matrices <- lapply(flatModel@matrices, applyDeps, pairs)
	flatModel@algebras <- lapply(flatModel@algebras, applyDeps, pairs)
	flatModel@objectives <- lapply(flatModel@objectives, applyDeps, pairs)
#	for (i in 1:length(edges)) {
#		source <- edgeNames[[i]]
#		sinks <- edges[[i]]
#		if (length(sinks) > 0) {
#			for (sink in sinks) {
#				flatModel[[sink]]@dependencies[[source]] <- 1L
#			}
#		}
#	}
	return(flatModel)
}

