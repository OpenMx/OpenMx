setClass(Class = "MxPathModel",
	representation = representation(
		paths = "data.frame",
		latentVars = "character",
		manifestVars = "character",
		matrices = "list"))

setMethod("initialize", "MxPathModel",
	function(.Object, paths = list(), latentVars = character(),
		manifestVars = character(), matrices = list()) {
		if (length(paths) > 0) {
			.Object <- mxAddPath(.Object, paths)
		}
		.Object@latentVars <- latentVars
		.Object@manifestVars <- manifestVars
		.Object@matrices <- matrices
		return(.Object)
	}
)

setGeneric("mxAddPath", function(.Object, paths) {
	return(standardGeneric("mxAddPath")) } )

setGeneric("mxRemovePath", function(.Object, paths) {
	return(standardGeneric("mxRemovePath")) } )	
	
froms <- function(lst) {
  retval <- lapply(lst, function(x) { return(x$from) } )
  return(retval)
}	

tos <- function(lst) {
  retval <- lapply(lst, function(x) { return(x$to) } )
  return(retval)
}	
	
setMethod("mxAddPath", "MxPathModel", 
	function(.Object, paths) {
		if (length(paths) < 1) {
			return(.Object)
		}
		if (isMxPath(paths)) {
			paths <- list(paths)
		}
		if (!all(sapply(paths, isMxPath))) {
			stop("Second argument is neither an MxPath nor a list of MxPaths")		
		}
		if (any(is.na(froms(paths))) || any(is.na(tos(paths)))) {		
			stop("The \'from\' field or the \'to\' field contains an NA")
		}
		for(i in 1:length(paths)) {
			.Object <- mxAddSinglePath(.Object, paths[[i]])
		}
		return(.Object)
	}
)

setMethod("mxRemovePath", "MxPathModel", 
	function(.Object, paths) {
		if (length(paths) < 1) {
			return(.Object)
		}
		if (isMxPath(paths)) {
			paths <- list(paths)
		}
		if (!all(sapply(paths, isMxPath))) {
			stop("Second argument is neither an MxPath nor a list of MxPaths")		
		}
		if (any(is.na(froms(paths))) || any(is.na(tos(paths)))) {		
			stop("The \'from\' field or the \'to\' field contains an NA")
		}		
		for(i in 1:length(paths)) {
			.Object <- mxRemoveSinglePath(.Object, paths[[i]])
		}
		return(.Object)
	}
)


mxAddSinglePath <- function(.Object, path) {
	if (nrow(.Object@paths) > 0) {
		fromExists <- (.Object@paths['from'] == path[['from']])
		toExists <- (.Object@paths['to'] == path[['to']])
		replace <- any(fromExists & toExists, na.rm=TRUE)
		morfExists <- (.Object@paths['from'] == path[['to']])
		otExists <- (.Object@paths['to'] == path[['from']])
		oppositeExists <- any(morfExists & otExists, na.rm=TRUE)
		if (oppositeExists) {
			newArrow <- !is.null(path[['arrows']])
			oldArrow <- !is.null(.Object@paths[morfExists & otExists,'arrows'])
			if (newArrow) {
				if (oldArrow && .Object@paths[morfExists & otExists,'arrows'] == 2) {
					fromTemp <- as.vector(.Object@paths[morfExists & otExists,'from'])
					toTemp <- as.vector(.Object@paths[morfExists & otExists,'to'])
					fUnique <- lapply(.Object@paths['from'], paste, collapse='')[[1]]
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

mxRemoveSinglePath <- function(.Object, path) {
	if (nrow(.Object@paths) > 0) {
		.Object@paths <- subset(.Object@paths, to != path[['to']] | from != path[['from']])
		if (nrow(.Object@paths) > 0) {		
			morfExists <- (.Object@paths['from'] == path[['to']])
			otExists <- (.Object@paths['to'] == path[['from']])
			oppositeExists <- any(morfExists & otExists, na.rm=TRUE)
			if (oppositeExists) {
				check1 <- !is.null(path[['arrows']]) && path[['arrows']] == 2
				check2 <- !is.null(.Object@paths[morfExists & otExists,'arrows']) &&
							.Object@paths[morfExists & otExists,'arrows'] == 2
				if (check1 || check2) {
					.Object@paths <- subset(.Object@paths, to != path[['from']] | from != path[['to']])
				}
			}
		}		
	}		
	return(.Object)
}

