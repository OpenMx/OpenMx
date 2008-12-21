setClass(Class = "MxModel",
	representation = representation(
		matrices = "list",
		algebras = "list",
		paths = "data.frame",
		latentVars = "character",
		manifestVars = "character",
		data = "data.frame"))
		
setMethod("initialize", "MxModel",
	function(.Object, paths = list(), latentVars = character(),
		manifestVars = character(), matrices = list(), 
		algebras = list(), data = data.frame()) {
		if (length(paths) > 0) {
			.Object <- mxAddPath(.Object, paths)
		}
		.Object@latentVars <- latentVars
		.Object@manifestVars <- manifestVars
		.Object@matrices <- matrices
		.Object@algebras <- algebras
		.Object@data <- data
		return(.Object)
	}
)

setMethod("[[", "MxModel",
	function(x, i, j, ..., drop = FALSE) {
		first <- x@matrices[[i]]
		second <- x@algebras[[i]]
		if (is.null(first)) {
			return(second)
		} else {
			return(first)
		}	
	}
)

setReplaceMethod("[[", "MxModel",
	function(x, i, j, value) {
		if (is(value,"MxMatrix")) {
			if (!is.null(x@algebras[[i]])) {
				stop(paste(i, "is already an MxAlgebra object"))
			}
			value@name <- i
			x@matrices[[i]] <- value
		} else if (is(value,"MxAlgebra")) {
			if (!is.null(x@matrices[[i]])) {
				stop(paste(i, "is already an MxMatrix object"))
			}
			value@name <- i
			x@algebras[[i]] <- value		
		} else {
			stop(paste("Unknown type of value", value))
		}
		return(x)
	}
)

setGeneric("omxAddEntries", function(.Object, entries) {
	return(standardGeneric("omxAddEntries")) } )

setGeneric("omxRemoveEntries", function(.Object, entries) {
	return(standardGeneric("omxRemoveEntries")) } )	
	
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
		model <- new("MxModel")	
	}
	lst <- list(...)
	if(class(model)[[1]] != "MxModel") {
		stop("First argument is not an MxModel object")
	}
	if(remove == TRUE) {
		model <- omxRemoveEntries(model, mappendHelper(lst, list()))
		if ( !is.null(manifestVars) ) {
			model@manifestVars <- setdiff(model@manifestVars, manifestVars)
		}
		if ( !is.null(latentVars) ) {
			model@latentVars <- setdiff(model@latentVars, latentVars)
		}				
	} else {
		model <- omxAddEntries(model, mappendHelper(lst, list()))
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

filterEntries <- function(entries, paths, matrices, algebras) {
	if (length(entries) == 0) {
		return(list(paths, matrices, algebras))
	}
	head <- entries[[1]]
	if (is(head, "MxMatrix")) {
		matrices <- append(matrices, head)
	} else if (is(head, "MxAlgebra")) {
		algebras <- append(algebras, head)
	} else if (omxIsPath(head)) {
		paths <- append(paths, head)
	} else {
		stop(paste("Unkown object:", head))
	}
	return(filterEntries(entries[-1], paths, matrices, algebras))
}
	
setMethod("omxAddEntries", "MxModel", 
	function(.Object, entries) {
		if (length(entries) < 1) {
			return(.Object)
		}
		threeTuple <- filterEntries(entries, list(), list(), list())
		paths <- threeTuple[[1]]
		matrices <- threeTuple[[2]]
		algebras <- threeTuple[[3]]
		if (any(is.na(froms(paths))) || any(is.na(tos(paths)))) {
			stop("The \'from\' field or the \'to\' field contains an NA")
		}
		if (length(paths) > 0) for(i in 1:length(paths)) {
			.Object <- omxAddSinglePath(.Object, paths[[i]])
		}
		if (length(matrices) > 0) for(i in 1:length(matrices)) {
			.Object <- omxAddSingleMatrix(.Object, matrices[[i]])
		}
		if (length(algebras) > 0) for(i in 1:length(algebras)) {
			.Object <- omxAddSingleAlgebra(.Object, algebras[[i]])
		}
		return(.Object)
	}
)

setMethod("omxRemoveEntries", "MxModel", 
	function(.Object, entries) {
		if (length(entries) < 1) {
			return(.Object)
		}
		threeTuple <- filterEntries(entries, list(), list(), list())
		paths <- threeTuple[[1]]
		matrices <- threeTuple[[2]]
		algebras <- threeTuple[[3]]		
		if (any(is.na(froms(paths))) || any(is.na(tos(paths)))) {
			stop("The \'from\' field or the \'to\' field contains an NA")
		}		
		if (length(paths) > 0) for(i in 1:length(paths)) {
			.Object <- omxRemoveSinglePath(.Object, paths[[i]])
		}
		if (length(matrices) > 0) for(i in 1:length(matrices)) {
			.Object <- omxRemoveSingleMatrix(.Object, matrices[[i]])
		}
		if (length(algebras) > 0) for(i in 1:length(algebras)) {
			.Object <- omxRemoveSingleAlgebra(.Object, algebras[[i]])
		}
		return(.Object)
	}
)

omxAddSingleMatrix <- function(.Object, matrix) {
	.Object[[matrix@name]] <- matrix
	return(.Object)
}

omxAddSingleAlgebra <- function(.Object, algebra) {
	.Object[[algebra@name]] <- algebra
	return(.Object)
}

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

omxRemoveSingleMatrix <- function(.Object, matrix) {
	.Object[[matrix@name]] <- NULL
	return(.Object)
}

omxRemoveSingleAlgebra <- function(.Object, algebra) {
	.Object[[algebra@name]] <- NULL
	return(.Object)
}

mxModel <- function(model = NULL, ..., remove = FALSE, manifestVars = NULL, latentVars = NULL) {
	omxModel(model, ..., remove = remove, manifestVars = manifestVars, latentVars = latentVars)
}

