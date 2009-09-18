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

cycleDetection <- function(flatModel) {
	dependencies <- new("MxDirectedGraph")
	objective <- flatModel@objective
	if (!is.null(objective) && !is.null(objective@dependencies)) {
		dependencies <- addObjectiveDetection(objective, dependencies)
	}
	if (length(flatModel@algebras) > 0) {
		for(i in 1:length(flatModel@algebras)) {
			dependencies <- addAlgebraDetection(flatModel@algebras[[i]], dependencies)
		}
		containsCycle(dependencies, flatModel)
	}
}

containsCycle <- function(graph, flatModel) {
	nodes <- graph@nodes
	if (length(nodes) > 0) {
		colors <- character()
		backedges <- character()
		colors[nodes] <- 'white'
		info <- list(colors, backedges)
		for(i in 1:length(nodes)) {
			node <- nodes[[i]]
			if (colors[[node]] == 'white') {
				info <- cycleVisitor(graph, node, info, flatModel@name)
				colors <- info[[1]]
			}
		}
	}
}

cycleVisitor <- function(graph, vertex, info, modelname) {
	colors <- info[[1]]
	backedges <- info[[2]]
	colors[[vertex]] <- 'grey'
	edges <- graph@edges[[vertex]]
	if (!is.null(edges) && length(edges) > 0) {
		for(i in 1:length(edges)) {
			destination <- edges[[i]]
			backedges[[destination]] <- vertex
			if (colors[[destination]] == 'grey') {
				reportCycle(backedges, destination, modelname)
			} else if (colors[[destination]] == 'white') {
				info <- list(colors, backedges)
				info <- cycleVisitor(graph, destination, info, modelname)
				colors <- info[[1]]
				backedges <- info[[2]]
			}
		}
	}
	colors[[vertex]] <- 'black'
	info <- list(colors, backedges)
	return(info)
}

reportCycle <- function(backedges, destination, modelname) {
	target <- backedges[[destination]]
	cycle <- union(destination, target)
	while(target != destination) {
		target <- backedges[[target]]
		cycle <- union(cycle, target)
	}
	cycle <- sapply(cycle, simplifyName, modelname)
	stop(paste("A cycle has been detected",
		"in model", omxQuotes(modelname),
		"involving the following elements:",
		omxQuotes(cycle)), call. = FALSE)
}

addObjectiveDetection <- function(objective, dependencies) {
	sources <- objective@dependencies
	sources <- sapply(sources, function(x) { slot(objective, x) })
	sources <- unlist(sources)
	sources <- sources[!is.na(sources)]
	sink <- objective@name
	if (length(sources) > 0) {
		dependencies <- omxAddDependency(sources, sink, dependencies)
	}
	return(dependencies)
}

addAlgebraDetection <- function(algebra, dependencies) {
	sink <- algebra@name
	formula <- algebra@formula
	dependencies <- addFormulaDetection(formula, sink, dependencies)
	return(dependencies)
}

omxAddDependency <- function(source, sink, dependencies) {
	dependencies <- addNode(source, dependencies)
	dependencies <- addNode(sink, dependencies)
	dependencies <- addEdge(source, sink, dependencies)
}

addFormulaDetection <- function(formula, sink, dependencies) {
	if (length(formula) == 1) {
		dependencies <- omxAddDependency(as.character(formula), sink, dependencies)
	} else {
		for (i in 2:length(formula)) {
			dependencies <- addFormulaDetection(formula[[i]], sink, dependencies)
		}
	}
	return(dependencies)
}

