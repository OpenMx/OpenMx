setClass(Class = "MxPathModel",
	representation = representation(
		paths = "data.frame",
		latentVars = "list",
		manifestVars = "list",
		matrices = "list"))

setMethod("initialize", "MxPathModel",
	function(.Object, paths = list(), latentVars = list(),
		manifestVars = list(), matrices = list()) {
#		.Object <- mxAddPath(.Object, paths)
		.Object@latentVars <- latentVars
		.Object@manifestVars <- manifestVars
		.Object@matrices <- matrices
		return(.Object)
	}
)

setGeneric("mxAddPath", function(.Object, paths) { return(standardGeneric("mxAddPath")) } )

setMethod("mxAddPath", "MxPathModel", 
	function(.Object, paths = list()) {
		if (isMxPath(paths)) {
			paths <- list(paths)
		}
		for(i in 1:length(paths)) {
			.Object <- mxAddSinglePath(.Object, paths[[i]])
		}
	}
)

mxAddSinglePath <- function(.Object, path) {
	if (nrow(.Object@paths) > 0) {
		fromExists <- (.Object@paths['from'] == path[['from']])
		toExists <- (.Object@paths['to'] == path[['to']])
		replace <- fromExists && toExists
		if (replace) {
			ids <- names(path)			
			for(i in 1:length(path)) {
				id <- ids[[i]]
				.Object@paths[fromExists & toExists,id] <- path[[id]]
			}
		} else {
			.Object@paths <- merge(.Object@paths, data.frame(path), all=TRUE)
		}
	} else {
		.Object@paths <- data.frame(path)
	}
	return(.Object)
}

