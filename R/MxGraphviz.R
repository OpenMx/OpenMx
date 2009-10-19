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


writeDotFile <- function(model, graph, dotFilename) {
	outputLines <- list()
	graphName <- paste('"', model@name, '"', sep='')
	outputLines <- c(outputLines, paste("digraph", graphName, "{"))
	outputLines <- c(outputLines, paste('\t', 'node [style=filled, fontname="Arial", fontsize=16];'))
	if (length(graph@manifestVars) > 0) {
		for(i in 1:length(graph@manifestVars)) {
			outputLines <- c(outputLines, paste('\t', graph@manifestVars[[i]], 
				'[shape=box, fillcolor="#a9fab1", height=0.5, width=0.5];'))
		}
	}
	if (length(graph@latentVars) > 0) {
		for(i in 1:length(graph@latentVars)) {
			outputLines <- c(outputLines, paste('\t', graph@latentVars[[i]], 
				'[shape=circle, fillcolor="#f4fd78"];'))		
		}
	}
	if (!is.null(model[['M']])) {
		outputLines <- c(outputLines, paste('\t', 'one', '[shape=triangle];'))
	}
	if (length(graph@paths) > 0) {
		for(i in 1:length(graph@paths)) {
			outputArrow <- list()
			path <- graph@paths[[i]]
			outputArrow <- c(outputArrow, paste('\t', path$from, "->", path$to))
			if (path$arrows == 1) {
				outputArrow <- c(outputArrow, "[dir=forward]")
			} else if (path$arrows == 2) {
			    if (path$from == path$to && path$from %in% graph@latentVars) {
					outputArrow <- c(outputArrow, "[dir=both, headport=n, tailport=n]")			
				} else if (path$from == path$to) {
					outputArrow <- c(outputArrow, "[dir=both, headport=s, tailport=s]")						
				} else {
					outputArrow <- c(outputArrow, "[dir=both]")					
				}
			}
			outputArrow <- c(outputArrow, ';')
			outputArrow <- paste(outputArrow, collapse='')
			outputLines <- c(outputLines, outputArrow)
		}
	}
	outputLines <- c(outputLines, '}\n')
	outputLines <- paste(outputLines, collapse='\n')
	cat(outputLines, file = dotFilename)
	return(invisible(outputLines))
}

omxGraphviz <- function(model, dotFilename = "") {
	if (missing(model) || !is(model, "MxModel")) {
		stop("The first argument is not an MxModel object")
	}
	if (!is(model, "MxRAMModel")) {
		stop(paste("The model", omxQuotes(model@name), 
			"is not a 'RAM' type model"))	
	}
	graph <- new("MxRAMGraph", model@manifestVars, model@latentVars)
	uniPaths <- matrixToPaths(model[['A']], 1)
	biPaths <- matrixToPaths(model[['S']], 2)
	meanPaths <- meansToPaths(model[['M']])
	graph@paths <- c(graph@paths, uniPaths, biPaths, meanPaths)
	writeDotFile(model, graph, dotFilename)
}

