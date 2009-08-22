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

graphvizAddMatrix <- function(values, rowFactors, colFactors, graph) {
	select <- (values != 0)
	if (length(select) > 0) {
		fromValues <- as.character(colFactors[select])
		toValues <- as.character(rowFactors[select])
		graph <- addEdge(fromValues, toValues, graph, 1)
	}
	return(graph)
}

omxGraphviz <- function(model) {
	loads <- require(Rgraphviz)
	if (!loads) return()
	if (missing(model) || !is(model, "MxModel")) {
		stop("The first argument is not an MxModel object")
	}
	if (!is(model, "MxRAMModel")) {
		stop(paste("The model", omxQuotes(model@name), "is not a 'RAM' type model"))	
	}
	graph <- new("graphNEL", edgemode = "directed")
	graph <- addNode(model@manifestVars, graph)
	graph <- addNode(model@latentVars, graph)

	aMatrix <- model[['A']]@values
	sMatrix <- model[['S']]@values
	rowFactors <- row(aMatrix, as.factor=TRUE)
	colFactors <- col(aMatrix, as.factor=TRUE)
	graph <- graphvizAddMatrix(aMatrix, rowFactors, colFactors, graph)
	graph <- graphvizAddMatrix(sMatrix, rowFactors, colFactors, graph)	
	
	return(graph)
}

