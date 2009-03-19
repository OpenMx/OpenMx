setClass(Class = "MxModel",
	representation = representation(
		name = "character",
		matrices = "list",
		algebras = "list",
		constraints = "list",
		bounds = "list",
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
		constraints = list(), bounds = list(),
		data = matrix(0, 0, 0), submodels = list(), 
		objective = NULL, independent = FALSE) {
		if (length(paths) > 0) {
			.Object <- mxAddPath(.Object, paths)
		}
		.Object@name <- name
		.Object@latentVars <- latentVars
		.Object@manifestVars <- manifestVars
		.Object@matrices <- matrices
		.Object@algebras <- algebras
		.Object@constraints <- constraints
		.Object@bounds <- bounds
		.Object@data <- data
		.Object@submodels <- submodels
		.Object@objective <- objective
		.Object@independent <- independent
		.Object@output <- list()
		return(.Object)
	}
)

omxExtractMethod <- function(model, index) {
	first <- model@matrices[[index]]
	second <- model@algebras[[index]]
	third <- model@submodels[[index]]
	fourth <- model@constraints[[index]]
	fifth <- model@bounds[[index]]
	if (!is.null(model@objective) && index == model@objective@name) {
		return(model@objective)
	} else if (!is.null(first)) {
		return(first)
	} else if (!is.null(second)) {
		return(second)
	} else if (!is.null(third)) {
		return(third)
	} else if (!is.null(fourth)) {
		return(fourth)
	} else {
		return(fifth)
	}	
}

omxReplaceMethod <- function(model, index, value) {
	current <- model[[index]]
	if (is.null(current) && is.null(value)) {
		return(model)
	}
	if(index == model@name) {
		stop(paste(omxQuotes(index), 
			"is already used as the name of the model"))
	}
	if(!is.null(current) && !is.null(value) && 
			!omxSameType(current, value)) {
		stop(paste("There already exists an object", 
				omxQuotes(index), 
				"in this model of different type"))
	}
	if(!is.null(value)) {
		value@name <- index
		test <- value		
	} else {
		test <- current
	}
	if (is(test,"MxMatrix")) {
		model@matrices[[index]] <- value
	} else if (is(test,"MxAlgebra")) {
		model@algebras[[index]] <- value
	} else if (is(test,"MxModel")) {
		model@submodels[[index]] <- value	
	} else if (is(test,"MxObjective")) {
		model@objective <- value
	} else if (is(test,"MxConstraint")) {
		model@constraints[[index]] <- value
	} else if (is(test,"MxBounds")) {
		model@bounds[[index]] <- value
	} else {
		stop("Unknown type of value", value)
	}
	return(model)
}

setMethod("[[", "MxModel",
	function(x, i, j, ..., drop = FALSE) {
		return(omxExtractMethod(x, i))
	}
)

setReplaceMethod("[[", "MxModel",
	function(x, i, j, value) {
		return(omxReplaceMethod(x, i, value))
	}
)

#setMethod("$", "MxModel",
#	function(x, name) {
#		return(omxExtractMethod(x, name))
#	}
#)
#
#setReplaceMethod("$", "MxModel",
#	function(x, name, value) {
#		return(omxReplaceMethod(x, name, value))
#	}
#)

omxSameType <- function(a, b) {
	return( (is(a, "MxModel") && is(b, "MxModel")) ||
			(is(a, "MxMatrix") && is(b, "MxMatrix")) ||
			(is(a, "MxAlgebra") && is(b, "MxAlgebra")) ||
			(is(a, "MxObjective") && is(b, "MxObjective")) ||
			(is(a, "MxConstraint") && is(b, "MxConstraint")) ||
			(is(a, "MxBounds") && is(b, "MxBounds")))
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

omxModel <- function(model = NA, ..., name = NA, manifestVars = NA,
	latentVars = NA, remove = FALSE, independent = NA) {
    first <- NULL
	if(typeof(model) != "S4" && is.na(model)) {
		model <- new("MxModel")	
	} else if (typeof(model) == "character") {
		model <- new("MxModel", name = model)
	} else if(!is(model, "MxModel")) {
 		first <- model
		model <- new("MxModel")
	}
	lst <- list(...)
	lst <- append(lst, first)
	if(class(model)[[1]] != "MxModel") {
		stop("First argument is not an MxModel object")
	}
	if(remove == TRUE) {
		model <- omxRemoveEntries(model, mappendHelper(lst, list()))
		if ( length(manifestVars) > 1 || !is.na(manifestVars) ) {
			model@manifestVars <- setdiff(model@manifestVars, manifestVars)
		}
		if ( length(latentVars) > 1 || !is.na(latentVars) ) {
			model@latentVars <- setdiff(model@latentVars, latentVars)
		}				
	} else {
		model <- omxAddEntries(model, mappendHelper(lst, list()))
		if ( length(manifestVars) > 1 || !is.na(manifestVars) ) {
			tmp <- append(model@manifestVars, manifestVars)
			model@manifestVars <- unique(tmp)
		}
		if (length(latentVars) > 1 || !is.na(latentVars)) {
			tmp <- append(model@latentVars, latentVars)
			model@latentVars <- unique(tmp)
		}		
	}
	if(!is.na(independent)) {
		model@independent <- independent
	}
	if(!is.na(name)) {
		model@name <- name
	}			
	return(model)
}

filterEntries <- function(entries, paths, namedEntities, 
	objectives, data) {
	if (length(entries) == 0) {
		return(list(paths, namedEntities,
			objectives, data))
	}
	head <- entries[[1]]
	pLength   <- length(paths)
	nLength   <- length(namedEntities)
	oLength   <- length(objectives)
	dLength   <- length(data)
	if (is.null(head)) {
    } else if (is(head, "MxObjective")) {
    	objectives[[oLength + 1]] <- head
    } else if (is(head, "MxData")) {
		data[[dLength + 1]] <- head
	} else if(omxIsPath(head)) {
		paths[[pLength + 1]] <- head
	} else if(isS4(head) && ("name" %in% slotNames(head))) {
		namedEntities[[nLength + 1]] <- head
	} else {
		stop(paste("Unknown object:", head))
	}
	return(filterEntries(entries[-1], paths, namedEntities, 
		objectives, data))
}
	
setMethod("omxAddEntries", "MxModel", 
	function(.Object, entries) {
		if (length(entries) < 1) {
			return(.Object)
		}
		tuple <- filterEntries(entries, list(), list(), list(), list())
		paths         <- tuple[[1]]
		namedEntities <- tuple[[2]]
		objectives    <- tuple[[3]]
		data          <- tuple[[4]]
		if (any(is.na(froms(paths))) || any(is.na(tos(paths)))) {
			stop("The \'from\' field or the \'to\' field contains an NA")
		}
		if (length(paths) > 0) for(i in 1:length(paths)) {
			.Object <- omxAddSinglePath(.Object, paths[[i]])
		}
		if (length(namedEntities) > 0) for(i in 1:length(namedEntities)) {
			.Object <- omxAddSingleNamedEntity(.Object, namedEntities[[i]])
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
		tuple <- filterEntries(entries, list(), list(), list(), list())
		paths         <- tuple[[1]]
		namedEntities <- tuple[[2]]
		objectives    <- tuple[[3]]
		data          <- tuple[[4]]
		if (any(is.na(froms(paths))) || any(is.na(tos(paths)))) {
			stop("The \'from\' field or the \'to\' field contains an NA")
		}		
		if (length(paths) > 0) for(i in 1:length(paths)) {
			.Object <- omxRemoveSinglePath(.Object, paths[[i]])
		}
		if (length(namedEntities) > 0) for(i in 1:length(namedEntities)) {
			.Object <- omxRemoveSingleNamedEntity(.Object, namedEntities[[i]])
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

omxAddSingleNamedEntity <- function(.Object, entity) {
	.Object[[entity@name]] <- entity
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

omxRemoveSingleNamedEntity <- function(.Object, entity) {
	.Object[[entity@name]] <- NULL
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
	cat("constraints :", omxQuotes(names(model@constraints)), '\n')
	cat("bounds :", omxQuotes(names(model@bounds)), '\n')
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

mxModel <- function(model = NA, ..., 
	manifestVars = NA, latentVars = NA, 
	remove = FALSE, independent = NA, name = NA) {
		omxModel(model, ..., name = name, 
		manifestVars = manifestVars, 
		latentVars = latentVars,
		remove = remove, 
		independent = independent)
}

