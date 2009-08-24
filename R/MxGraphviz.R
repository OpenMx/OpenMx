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

setClass(Class = "MxGraph",
	representation = representation(
		manifestVars = "character", 
		latentVars = "character", 
		paths = "list"))

setMethod("initialize", "MxGraph",
	function(.Object, manifestVars = character(), 
		latentVars = character(), paths = list()) {
		.Object@manifestVars <- manifestVars
		.Object@latentVars <- latentVars
		.Object@paths <- paths
		return(.Object)
	}
)	

writeDotFile <- function(model, graph, dotFilename) {
	dotFile <- file(dotFilename, "w")
	graphName <- sub(" ", "_", model@name, fixed=TRUE)
	cat("digraph", graphName, "{", "\n", file = dotFile)
	cat('\t', 'node [style=filled,fontname="Arial",fontsize=16]\n', file = dotFile)
	if (length(graph@manifestVars) > 0) {
		for(i in 1:length(graph@manifestVars)) {
			cat('\t', graph@manifestVars[[i]], 
				'[shape=box,fillcolor="#a9fab1",height=0.5,width=0.5];\n', 
				file = dotFile)
		}
	}
	if (length(graph@latentVars) > 0) {
		for(i in 1:length(graph@latentVars)) {
			cat('\t', graph@latentVars[[i]], 
				'[shape=circle,fillcolor="#f4fd78"];\n', file = dotFile)		
		}
	}
	if (length(graph@paths) > 0) {
		for(i in 1:length(graph@paths)) {
			path <- graph@paths[[i]]
			cat('\t', path$from, "->", path$to, file = dotFile)
			if (path$arrows == 1) {
				cat("[dir=forward]", file = dotFile)
			} else if (path$arrows == 2) {
			    if (path$from == path$to && path$from %in% graph@latentVars) {
					cat("[dir=both,headport=n,tailport=n]", file = dotFile)			
				} else if (path$from == path$to) {
					cat("[dir=both,headport=s,tailport=s]", file = dotFile)						
				} else {
					cat("[dir=both;]", file = dotFile)					
				}
			}
			cat(';\n', file = dotFile)
		}
	}
	cat("}", "\n", file = dotFile)
	close(dotFile)
}

omxGraphviz <- function(model, dotFilename) {
	if (missing(model) || !is(model, "MxModel")) {
		stop("The first argument is not an MxModel object")
	}
	if (!is(model, "MxRAMModel")) {
		stop(paste("The model", omxQuotes(model@name), 
			"is not a 'RAM' type model"))	
	}
	graph <- new("MxGraph", model@manifestVars, model@latentVars)
	uniPaths <- matrixToPaths(model[['A']], 1)
	biPaths <- matrixToPaths(model[['S']], 2)
	graph@paths <- c(graph@paths, uniPaths, biPaths)
	writeDotFile(model, graph, dotFilename)
}

