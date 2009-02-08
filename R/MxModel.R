setClass(Class = "MxModel",
	representation = representation(
		name = "character",
		matrices = "list",
		algebras = "list",
		paths = "data.frame",
		latentVars = "character",
		manifestVars = "character",
		data = "MxData",
		submodels = "list",
		objective = "MxObjective",
		independent = "logical",
		output = "list"
))

setMethod("initialize", "MxModel",
	function(.Object, name = omxUntitledName(), paths = list(), 
		latentVars = character(), manifestVars = character(), 
		matrices = list(), algebras = list(), 
		data = NULL, submodels = list(), 
		objective = NULL, independent = FALSE) {
		if (length(paths) > 0) {
			.Object <- mxAddPath(.Object, paths)
		}
		.Object@name <- name
		.Object@latentVars <- latentVars
		.Object@manifestVars <- manifestVars
		.Object@matrices <- matrices
		.Object@algebras <- algebras
		.Object@data <- data
		.Object@submodels <- submodels
		.Object@objective <- objective
		.Object@independent <- independent
		.Object@output <- list()
		return(.Object)
	}
)

setMethod("[[", "MxModel",
	function(x, i, j, ..., drop = FALSE) {
		first <- x@matrices[[i]]
		second <- x@algebras[[i]]
		third <- x@submodels[[i]]
		if (!is.null(x@objective) && i == x@objective@name) {
			return(x@objective)
		} else if (is.null(first) && is.null(second)) {
			return(third)
		} else if (is.null(first)) {
			return(second)
		} else {
			return(first)
		}
	}
)

setReplaceMethod("[[", "MxModel",
	function(x, i, j, value) {
		current <- x[[i]]
		if (is.null(current) && is.null(value)) {
			return(x)
		}
		if(i == x@name) {
			stop(paste(omxQuotes(i), 
				"is already used as the name of the model"))
		}
		if(!is.null(current) && !is.null(value) && 
				!omxSameType(current, value)) {
			stop(paste("There already exists an object", 
					omxQuotes(i), 
					"in this model of different type"))
		}
		if(!is.null(value)) {
			value@name <- i
			test <- value		
		} else {
			test <- current
		}
		if (is(test,"MxMatrix")) {
			x@matrices[[i]] <- value
		} else if (is(test,"MxAlgebra")) {
			x@algebras[[i]] <- value		
		} else if (is(test,"MxModel")) {
			x@submodels[[i]] <- value			
		} else if (is(test,"MxObjective")) {
			x@objective <- value
		} else {
			stop("Unknown type of value", value)
		}
		return(x)
	}
)

omxSameType <- function(a, b) {
	return( (is(a, "MxModel") && is(b, "MxModel")) ||
			(is(a, "MxMatrix") && is(b, "MxMatrix")) ||
			(is(a, "MxAlgebra") && is(b, "MxAlgebra")) ||
			(is(a, "MxObjective") && is(b, "MxObjective")))
}

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
		len <- length(result)
		addition <- lst[[1]]
		if (is.list(addition)) {
			result <- append(result, addition)
		} else {
			result[[len + 1]] <- addition
		}
		return(result)
	} else {
		len <- length(result)
		addition <- lst[[1]]
		if (is.list(addition)) {
			result <- append(result, addition)
		} else {
			result[[len + 1]] <- addition
		}
		return(mappendHelper(lst[2:length(lst)], result))
	}
}

omxModel <- function(model = NULL, ..., name = NULL, manifestVars = NULL,
	latentVars = NULL, remove = FALSE, independent = NULL) {
	if(is.null(model)) {
		model <- new("MxModel")	
	} else if (typeof(model) == "character") {
		model <- new("MxModel", name = model)
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
	if(!is.null(independent)) {
		model@independent <- independent
	}
	if(!is.null(name)) {
		model@name <- name
	}			
	return(model)
}

filterEntries <- function(entries, paths, matrices, algebras, models, objectives, data) {
	if (length(entries) == 0) {
		return(list(paths, matrices, algebras, models, objectives, data))
	}
	head <- entries[[1]]
	pLength   <- length(paths)
	matLength <- length(matrices)
	aLength   <- length(algebras)	
	modLength <- length(models)
	oLength   <- length(objectives)
	dLength   <- length(data)
	if (is.null(head)) {
	} else if (is(head, "MxMatrix")) {
		matrices[[matLength + 1]] <- head
	} else if (is(head, "MxAlgebra")) {
		algebras[[aLength + 1]] <- head
	} else if (omxIsPath(head)) {		
		paths[[pLength + 1]] <- head
    } else if (is(head, "MxModel")) {
    	models[[modLength + 1]] <- head
    } else if (is(head, "MxObjective")) {
    	objectives[[oLength + 1]] <- head
    } else if (is(head, "MxData")) {
		data[[dLength + 1]] <- head
	} else if(omxIsPath(head)) {
		paths[[pLength + 1]] <- head
	} else {
		stop(paste("Unknown object:", head))
	}
	return(filterEntries(entries[-1], paths, matrices, algebras, models, objectives, data))
}
	
setMethod("omxAddEntries", "MxModel", 
	function(.Object, entries) {
		if (length(entries) < 1) {
			return(.Object)
		}
		tuple <- filterEntries(entries, list(), list(), list(), list(), list(), list())
		paths      <- tuple[[1]]
		matrices   <- tuple[[2]]
		algebras   <- tuple[[3]]
		models     <- tuple[[4]]
		objectives <- tuple[[5]]
		data       <- tuple[[6]]
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
		if (length(models) > 0) for(i in 1:length(models)) {
			.Object <- omxAddSingleModel(.Object, models[[i]])
		}
		if (length(objectives) > 0) .Object <- omxAddObjectives(.Object, objectives)
		if (length(data) > 0) .Object <- omxAddData(.Object, data)
		return(.Object)
	}
)

setMethod("omxRemoveEntries", "MxModel", 
	function(.Object, entries) {
		if (length(entries) < 1) {
			return(.Object)
		}
		tuple <- filterEntries(entries, list(), list(), list(), list(), list(), list())
		paths      <- tuple[[1]]
		matrices   <- tuple[[2]]
		algebras   <- tuple[[3]]
		models     <- tuple[[4]]
		objectives <- tuple[[5]]
		data       <- tuple[[6]]
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
		if (length(models) > 0) for(i in 1:length(models)) {
			.Object <- omxRemoveSingleModel(.Object, models[[i]])
		}
		if (length(objectives) > 0) {
			.Object@objective <- NULL
		}
		if (length(data) > 0) {
			.Object@data <- NULL
		}
		return(.Object)
	}
)

omxAddSingleModel <- function(.Object, model) {
	.Object[[model@name]] <- model
	return(.Object)	
}

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

omxRemoveSingleModel <- function(.Object, model) {
	.Object[[model@name]] <- NULL
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

omxAddObjectives <- function(.Object, objectives) {
	if (length(objectives) > 1) {
		warning("Multiple objective functions were specified; the first one will be used")
	}
	objective <- objectives[[1]]
	.Object[[objective@name]] <- objective
	return(.Object)
}

omxAddData <- function(.Object, dataset) {
	if (length(dataset) > 1) {
		warning("Multiple datasets were specified; the first one will be used")
	}
	data <- dataset[[1]]
	.Object@data <- data
	return(.Object)
}

omxQuotes <- function(name) {
	listTerms <- sapply(name, function(x) {paste("'", x, "'", sep = '')} )
	return(paste(listTerms, collapse=', '))
}

omxDisplayModel <- function(model) {
	cat("MxModel", omxQuotes(model@name), '\n')
	cat("matrices :", omxQuotes(names(model@matrices)), '\n')
	cat("algebras :", omxQuotes(names(model@algebras)), '\n')
	if (length(model@paths) > 0) {
		cat("latentVars :", model@latentVars, '\n')
		cat("manifestVars :", model@manifestVars, '\n')
		cat("paths :", nrow(model@paths), "paths", '\n')
	}
	if (is.null(model@data)) {
		cat("data : NULL\n")
	} else {
		cat("data :", nrow(model@data), "x", ncol(model@data), '\n')
	}
	cat("submodels :", omxQuotes(names(model@submodels)), '\n')
	objective <- model@objective
	if (is.null(objective)) {
		objectiveType <- "NULL"
		objectiveName <- ""
	} else {
		objectiveType <- class(objective)[[1]]
		objectiveName <- omxQuotes(objective@name)
	}
	cat("objective :", objectiveType, objectiveName, '\n')
	cat("independent :", model@independent, '\n')
	cat("output :", length(model@output) > 0, '\n')
}

setMethod("print", "MxModel", function(x,...) { omxDisplayModel(x) })
setMethod("show", "MxModel", function(object) { omxDisplayModel(object) })

mxModel <- function(model = NULL, ..., name = NULL, manifestVars = NULL,
	latentVars = NULL, remove = FALSE, independent = NULL) {
	omxModel(model, ..., name = name, manifestVars = manifestVars, latentVars = latentVars,
		remove = remove, independent = independent)
}

