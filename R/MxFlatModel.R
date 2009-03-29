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


setClass(Class = "MxFlatModel",
	representation = representation(
		objectives = "list",
		datasets = "list"
	),
	contains = "MxModel")
	
setMethod("initialize", "MxFlatModel",
	function(.Object, model, objectives, datasets) {
		modelSlotNames <- slotNames(model)
		for(i in 1:length(modelSlotNames)) {
			name <- modelSlotNames[[i]]
			slot(.Object, name) <- slot(model, name)
		}
		.Object@objectives <- objectives
		.Object@datasets <- datasets
		return(.Object)
	}
)	

omxGenerateDefinitionNames <- function(datasets) {
	nameList <- lapply(datasets, 
		function(x) { dimnames(x@matrix)[[2]] })
	result <- list()
	if(length(nameList) > 0) {
		for(i in 1:length(nameList)) {
			colNames <- nameList[[i]]
			if(length(colNames) > 0) {
				for(j in 1:length(colNames)) {
					name <- colNames[[j]]
					result[[name]] <- c(i - 1, j - 1)
				}
			}	
		}
	}
	return(result)
}


setMethod("print", "MxFlatModel", function(x,...) {
	callNextMethod()
	cat("objectives : ")
	print(x@objectives)
	cat("datasets :", length(x@datasets), '\n') 
})

setMethod("show", "MxFlatModel", function(object) { 
	callNextMethod()
	cat("objectives : ")
	print(object@objectives)
	cat("datasets :", length(object@datasets), '\n') 
})
