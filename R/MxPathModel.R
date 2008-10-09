setClass(Class = "MxPathModel",
	representation = representation(
		paths = "data.frame",
		latentVars = "list",
		manifestVars = "list",
		matrices = "list"))

setMethod("initialize", "MxPathModel",
	function(.Object, paths = list(), latentVars = list(),
		manifestVars = list(), matrices = list()) {
		.Object <- mxAddPath(.Object, paths)
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
		replace <- any(fromExists & toExists)
		morfExists <- (.Object@paths['from'] == path[['to']])
		otExists <- (.Object@paths['to'] == path[['from']])
		oppositeExists <- any(morfExists & otExists)
		if (oppositeExists) {
			newArrow <- !is.null(path[['arrows']])
			oldArrow <- !is.null(.Object@paths[morfExists & otExists,'arrows'])
			if (newArrow && oldArrow) {
				if (.Object@paths[morfExists & otExists,'arrows'] == 2) {
					fromTemp <- as.vector(.Object@paths[morfExists & otExists,'from'])
					toTemp <- as.vector(.Object@paths[morfExists & otExists,'to'])
					fUnique <- lapply(.Object@paths['from'], paste)[[1]]
					.Object@paths[morfExists & otExists, 'from'] <- fUnique
					.Object@paths[.Object@paths['from'] == fUnique, 'to'] <- fromTemp
					.Object@paths[.Object@paths['from'] == fUnique, 'from'] <- toTemp
					fromExists <- (.Object@paths['from'] == path[['from']])
					toExists <- (.Object@paths['to'] == path[['to']])
					replace <- TRUE
				} else if (path[['arrows']] == 2) {
					tmp <- path[['from']]
					path[['from']] <- path[['to']]
					path[['to']] <- tmp
					fromExists <- (.Object@paths['from'] == path[['from']])
					toExists <- (.Object@paths['to'] == path[['to']])
					replace <- TRUE
				}
			}
		}
		if (replace) {
			ids <- names(path)			
			for(i in 1:length(path)) {
				id <- ids[[i]]
				.Object@paths[fromExists & toExists,id] <- path[[id]]
			}
		} else {
			.Object@paths <- merge(.Object@paths, 
				data.frame(path, stringsAsFactors = FALSE), all=TRUE)
		}
	} else {
		.Object@paths <- data.frame(path, stringsAsFactors = FALSE)
	}
	return(.Object)
}

