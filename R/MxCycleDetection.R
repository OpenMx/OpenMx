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

cycleDetection <- function(flatModel) {
	dependencies <- new("MxDirectedGraph")
	if (length(flatModel@fitfunctions) > 0) {
		for(i in 1:length(flatModel@fitfunctions)) {
			dependencies <- addFitFunctionDetection(flatModel@fitfunctions[[i]], flatModel, dependencies)
		}	
	}
	if (length(flatModel@expectations) > 0) {
		for(i in 1:length(flatModel@expectations)) {
			dependencies <- addExpectationDetection(flatModel@expectations[[i]], dependencies)
		}	
	}
	if (length(flatModel@algebras) > 0) {
		for(i in 1:length(flatModel@algebras)) {
			dependencies <- addAlgebraDetection(flatModel@algebras[[i]], dependencies)
		}
	}
	if (length(flatModel@matrices) > 0) {
		for(i in 1:length(flatModel@matrices)) {
			dependencies <- addMatrixDetection(flatModel@matrices[[i]], dependencies)
		}
	}
	containsCycle(dependencies, flatModel)
	return(dependencies)
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
	report <- cycle[!sapply(cycle, hasSquareBrackets)]
	stop(paste("A cycle has been detected",
		"in model", omxQuotes(modelname),
		"involving the following elements:",
		omxQuotes(report)), call. = FALSE)
}

addFitFunctionDetection <- function(fitfunction, flatModel, dependencies) {
	dependencies <- genericFitDependencies(fitfunction, flatModel, dependencies)
	return(dependencies)
}

addExpectationDetection <- function(expectation, dependencies) {
	dependencies <- genericExpDependencies(expectation, dependencies)
	return(dependencies)
}

addMatrixDetection <- function(matrix, dependencies) {
	labels <- matrix@labels
	select <- matrix@.squareBrackets
	subs <- labels[select]
	if (length(subs) > 0) {
		for(i in 1:length(subs)) {
			components <- splitSubstitution(subs[[i]])
			name <- components[[1]]
			dependencies <- imxAddDependency(name, matrix@name, dependencies)
		}
	}
	return(dependencies)
}

addAlgebraDetection <- function(algebra, dependencies) {
	if (algebra@fixed) return(dependencies)
	sink <- algebra@name
	formula <- algebra@formula
	dependencies <- addFormulaDetection(formula, sink, dependencies)
	return(dependencies)
}

addFormulaDetection <- function(formula, sink, dependencies) {
	if (length(formula) == 1) {
		dependencies <- imxAddDependency(as.character(formula), sink, dependencies)
	} else {
		for (i in 2:length(formula)) {
			dependencies <- addFormulaDetection(formula[[i]], sink, dependencies)
		}
	}
	return(dependencies)
}

##' Add a dependency
##'
##' The dependency tracking system ensures that algebra and
##' fitfunctions are not recomputed if their inputs have not changed.
##' Dependency information is computed prior to handing the model off
##' to the optimizer to reduce overhead during optimization.
##'
##' Each free parameter keeps track of all the objects that store that
##' free parameter and the transitive closure of all algebras and fit
##' functions that depend on that free parameter.  Similarly, each
##' definition variable keeps track of all the objects that store that
##' free parameter and the transitive closure of all the algebras and
##' fit functions that depend on that free parameter. At each
##' iteration of the optimization, when the free parameter values are
##' updated, all of the dependencies of that free parameter are marked
##' as dirty (see \code{omxFitFunction.repopulateFun}). After an
##' algebra or fit function is computed, \code{omxMarkClean()} is
##' called to to indicate that the algebra or fit function is updated.
##' Similarly, when definition variables are populated in FIML, all of
##' the dependencies of the definition variables are marked as dirty.
##' Particularly for FIML, the fact that non-definition-variable
##' dependencies remain clean is a big performance gain.
##'
##' @param source a character vector of the names of the computation sources (inputs)
##' @param sink the name of the computation sink (output)
##' @param dependencies the dependency graph

imxAddDependency <- function(source, sink, dependencies) {
	if (length(source) == 0) {
		warning("imxAddDependency called with no sources (ignored)")
		return(dependencies)
	}
	dependencies <- addNode(source, dependencies)
	dependencies <- addNode(sink, dependencies)
	dependencies <- addEdge(source, sink, dependencies)
}

