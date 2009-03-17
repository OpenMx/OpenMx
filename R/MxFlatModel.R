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
		function(x) { dimnames(x)[[2]] })
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