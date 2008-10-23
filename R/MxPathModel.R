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

setGeneric("omxAddPath", function(.Object, paths) {
	return(standardGeneric("omxAddPath")) } )

setGeneric("omxRemovePath", function(.Object, paths) {
	return(standardGeneric("omxRemovePath")) } )	
	
froms <- function(lst) {
  retval <- lapply(lst, function(x) { return(x$from) } )
  return(retval)
}	

tos <- function(lst) {
  retval <- lapply(lst, function(x) { return(x$to) } )
  return(retval)
}	

omxMappend <- function(...) {
    args <- list(...)
	return(mappendHelper(args, list()))
}

mappendHelper <- function(lst, result) {
	if (length(lst) == 0) {
		return(result)
	} else if (length(lst) == 1) {
		return(append(result,lst[[1]]))
	} else {
		return(mappendHelper(lst[2:length(lst)], append(result, lst[[1]])))
	}
}

omxModel <- function(model=NULL, ..., remove=FALSE, manifestVars=NULL, latentVars=NULL) {
	if(is.null(model)) {
		model <- new("MxPathModel")	
	}
	lst <- list(...)
	if(class(model)[[1]] != "MxPathModel") {
		stop("First argument is not an MxPathModel object")
	}
	if(remove == TRUE) {
		model <- omxRemovePath(model, mappendHelper(lst, list()))
		if ( !is.null(manifestVars) ) {
			model@manifestVars <- setdiff(model@manifestVars, manifestVars)
		}
		if ( !is.null(latentVars) ) {
			model@latentVars <- setdiff(model@latentVars, latentVars)
		}				
	} else {
		model <- omxAddPath(model, mappendHelper(lst, list()))
		if ( !is.null(manifestVars) ) {
			tmp <- append(model@manifestVars, manifestVars)
			model@manifestVars <- unique(tmp)
		}
		if ( !is.null(latentVars) ) {
			tmp <- append(model@latentVars, latentVars)
			model@latentVars <- unique(tmp)
		}		
	}
	return(model)
}
	
setMethod("omxAddPath", "MxPathModel", 
	function(.Object, paths) {
		if (length(paths) < 1) {
			return(.Object)
		}
		if (!all(sapply(paths, omxIsPath))) {
			stop("Second argument is not a list of MxPaths")
		}
		if (any(is.na(froms(paths))) || any(is.na(tos(paths)))) {
			stop("The \'from\' field or the \'to\' field contains an NA")
		}
		for(i in 1:length(paths)) {
			.Object <- omxAddSinglePath(.Object, paths[[i]])
		}
		return(.Object)
	}
)

setMethod("omxRemovePath", "MxPathModel", 
	function(.Object, paths) {
		if (length(paths) < 1) {
			return(.Object)
		}
		if (!all(sapply(paths, omxIsPath))) {
			stop("Second argument is not a list of MxPaths")
		}
		if (any(is.na(froms(paths))) || any(is.na(tos(paths)))) {
			stop("The \'from\' field or the \'to\' field contains an NA")
		}		
		for(i in 1:length(paths)) {
			.Object <- omxRemoveSinglePath(.Object, paths[[i]])
		}
		return(.Object)
	}
)


omxAddSinglePath <- function(.Object, path) {
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
			} else if (newArrow && path[['arrows']] == 2) {
					tmp <- path[['from']]
					path[['from']] <- path[['to']]
					path[['to']] <- tmp
					fromExists <- (.Object@paths['from'] == path[['from']])
					toExists <- (.Object@paths['to'] == path[['to']])
					replace <- TRUE
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
	fromExists <- (.Object@paths['from'] == path[['from']])
	toExists <- (.Object@paths['to'] == path[['to']])
	field <- .Object@paths[fromExists & toExists, 'arrows']
	if (!is.null(field) && !is.na(field)  
			&& (field == 2) 
			&& (path[['from']] > path[['to']])) {
		fromTemp <- as.vector(.Object@paths[morfExists & otExists,'from'])
		toTemp <- as.vector(.Object@paths[morfExists & otExists,'to'])
		fUnique <- lapply(.Object@paths['from'], paste, collapse='')[[1]]
		.Object@paths[morfExists & otExists, 'from'] <- fUnique
		.Object@paths[.Object@paths['from'] == fUnique, 'to'] <- fromTemp
		.Object@paths[.Object@paths['from'] == fUnique, 'from'] <- toTemp
	}	
	return(.Object)
}

omxRemoveSinglePath <- function(.Object, path) {
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
					.Object@paths <- subset(.Object@paths, 
						to != path[['from']] | from != path[['to']])
				}
			}
		}		
	}		
	return(.Object)
}

mxModel <- function(model=NULL, ..., remove=FALSE, manifestVars=NULL, latentVars=NULL) {	omxModel(model, ..., remove = remove, manifestVars = manifestVars, latentVars = latentVars)
}
