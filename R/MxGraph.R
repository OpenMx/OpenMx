#
#   Copyright 2007-2016 The OpenMx Project
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

##' @title MxRAMGraph
##' @name MxRAMGraph-class
##'
##' @description
##' This is an internal class and should not be used directly.
##' It is a class for RAM directed graphs.
##'
##' @aliases
##' MxRAMGraph
##' @rdname MxRAMGraph-class
setClass(Class = "MxRAMGraph",
	representation = representation(
		manifestVars = "character", 
		latentVars = "character", 
		paths = "list"))

setMethod("initialize", "MxRAMGraph",
	function(.Object, manifestVars = character(), 
		latentVars = character(), paths = list()) {
		.Object@manifestVars <- manifestVars
		.Object@latentVars <- latentVars
		.Object@paths <- paths
		return(.Object)
	}
)

##' @title MxDirectedGraph
##' @name MxDirectedGraph-class
##'
##' @description
##' This is an internal class and should not be used directly.
##' It is a class for directed graphs.
##'
##' @aliases
##' MxDirectedGraph
##' @rdname MxDirectedGraph-class
setClass(Class = "MxDirectedGraph",
	representation = representation(
		nodes = "character", 
		edges = "list"))

setMethod("initialize", "MxDirectedGraph",
	function(.Object, nodes = character(), 
		edges = list()) {
		.Object@nodes <- nodes
		.Object@edges <- edges
		return(.Object)
	}
)

setGeneric("addNode", function(nodes, graph) {
	return(standardGeneric("addNode")) } )

setGeneric("addEdge", function(source, sink, graph) {
	return(standardGeneric("addEdge")) } )

setMethod("addNode", c("character", "MxDirectedGraph"), 
	function(nodes, graph) { 
	if (any(is.na(nodes))) {
		stop("NA is not a valid node name")
	}
	graph@nodes <- union(graph@nodes, nodes)
	return(graph)
})

setMethod("addEdge", c("character", "character", "MxDirectedGraph"), 
	function(source, sink, graph) { 
	if (any(is.na(source)) || any(is.na(sink))) {
		stop("NA is not a valid node name")
	}
	if (length(source) == 0) {
		stop("You must specify at least one source")
	}
	if (length(source) == 0) {
		stop("You must specify at least one sink")
	}
	missingsources <- setdiff(source, graph@nodes)
	missingsinks <- setdiff(sink, graph@nodes)
	missing <- union(missingsources, missingsinks)
	if (length(missing) > 0) {
		stop(paste("The following sources/sinks",
			"are not nodes in the graph:",
			omxQuotes(missing)))
	}
	pairs <- mapply(function(x,y) { c(x,y) }, source, sink, SIMPLIFY = FALSE)
	edges <- graph@edges	
	for (i in 1:length(pairs)) {
		pair <- pairs[[i]]
		asource <- pair[[1]]
		asink <- pair[[2]]
		target <- edges[[asource]]
		if (is.null(target)) {
			edges[[asource]] <- asink
		} else {
			target <- union(target, asink)
			edges[[asource]] <- target
		}
	}
	graph@edges <- edges
	return(graph)
})
