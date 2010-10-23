#
#   Copyright 2007-2010 The OpenMx Project
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

setClass(Class = "MxRAMModel",
	representation = representation(),
	contains = "MxModel")

omxModelTypes[['RAM']] <- "MxRAMModel"

# Define generic functions

setMethod("omxTypeName", "MxRAMModel", 
	function(model) { "RAM" }
)

setMethod("omxInitModel", "MxRAMModel", 
	function(model) {
		if (is.null(model[['A']])) {
			model[['A']] <- createMatrixA(model)
		}
		if (is.null(model[['S']])) {
			model[['S']] <- createMatrixS(model)
		}
		if (is.null(model[['objective']])) {
			model[['objective']] <- mxRAMObjective('A', 'S', 'F')
		}
		model[['F']] <- createMatrixF(model)
		return(model)
	}
)

setMethod("omxModelBuilder", "MxRAMModel", 
	function(model, lst, name, 
		manifestVars, latentVars, remove, independent) {
		model <- variablesArgumentRAM(model, manifestVars, latentVars, remove)
		model <- listArgumentRAM(model, lst, remove)
		notPathOrData <- getNotPathsOrData(lst)
		callNextMethod(model, notPathOrData, name, character(), 
			character(), remove, independent)
	}
)

setMethod("omxVerifyModel", "MxRAMModel",
    function(model) {
        if ((length(model$A) == 0) ||
             (length(model$S) == 0) ||
             (length(model$F) == 0)) {
			msg <- paste("The RAM model", omxQuotes(model@name),
                "does not contain any paths.")
			stop(msg, call. = FALSE)
        }
        objective <- model$objective
        if (!is.null(objective) && !is(objective, "MxRAMObjective")) {
			msg <- paste("The RAM model", omxQuotes(model@name),
                "does not contain a RAM objective function.")
			stop(msg, call. = FALSE)        	
        }
        if (!is.null(model@data) && model@data@type == "raw" &&
        	is.null(model$M)) {
			msg <- paste("The RAM model", omxQuotes(model@name),
                "contains raw data but has not specified any means paths.")
			stop(msg, call. = FALSE)
        }
        if (!is.null(model@data) && !single.na(model@data@means) &&
        	is.null(model$M)) {
			msg <- paste("The RAM model", omxQuotes(model@name),
                "contains an observed means vector",
                "but has not specified any means paths.")
			stop(msg, call. = FALSE)        	
		}
        if (length(model@submodels) > 0) {
        	return(all(sapply(model@submodels, omxVerifyModel)))
        }        
        return(TRUE)
    }
)


setReplaceMethod("[[", "MxRAMModel",
	function(x, i, j, value) {
		return(replaceMethodRAM(x, i, value))
	}
)

setReplaceMethod("$", "MxRAMModel",
	function(x, name, value) {
		return(replaceMethodRAM(x, name, value))
	}
)

# Helper functions used by the generic functions

variablesArgumentRAM <- function(model, manifestVars, latentVars, remove) {
	if (single.na(manifestVars)) {
		manifestVars <- character()
	}
	if (single.na(latentVars)) {
		latentVars <- character()
	}
	if (remove == TRUE) {
		if (length(latentVars) + length(manifestVars) > 0) {
			model <- removeVariablesRAM(model, latentVars, manifestVars)
		}
	} else if (length(manifestVars) + length(latentVars) > 0) {
		latentVars <- as.character(latentVars)
		manifestVars <- as.character(manifestVars)
		checkVariables(model, latentVars, manifestVars)
		model <- addVariablesRAM(model, latentVars, manifestVars)
	}
	return(model)
}

removeVariablesRAM <- function(model, latent, manifest) {
	missingLatent <- setdiff(latent, model@latentVars)
	missingManifest <- setdiff(manifest, model@manifestVars)
	if (length(missingLatent) > 0) {
		stop(paste("The latent variable(s)", omxQuotes(missingLatent),
			"are not present in the model.",
			"They cannot be deleted"), call. = FALSE)
	} else if (length(missingManifest) > 0) {
		stop(paste("The manifest variable(s)", omxQuotes(missingManifest),
			"are not present in the model.",
			"They cannot be deleted"), call. = FALSE)
	} else if (length(unique(latent)) != length(latent)) {
		stop("The latent variables list contains duplicate elements",
			call. = FALSE)
	} else if (length(unique(manifest)) != length(manifest)) {
		stop("The manifest variables list contains duplicate elements",
			call. = FALSE)
	}
	model@latentVars <- setdiff(model@latentVars, latent)
	model@manifestVars <- setdiff(model@manifestVars, manifest)
	A <- model[['A']]
	S <- model[['S']]
	A <- removeVariablesAS(A, latent)
	A <- removeVariablesAS(A, manifest)
	S <- removeVariablesAS(S, latent)
	S <- removeVariablesAS(S, manifest)
	model[['F']] <- createMatrixF(model)
	model[['A']] <- A
	model[['S']] <- S
	return(model)
}

addVariablesRAM <- function(model, latent, manifest) {
	newLatent   <- setdiff(latent, model@latentVars)
	newManifest <- setdiff(manifest, model@manifestVars)
	model@latentVars <- c(model@latentVars, newLatent)
	model@manifestVars <- c(model@manifestVars, newManifest)
	A <- model[['A']]
	S <- model[['S']]
	M <- model[['M']]
	if (is.null(A)) {
		A <- createMatrixA(model)
	} else {
		A <- addVariablesAS(A, model, newLatent, newManifest)
	}
	if (is.null(S)) {
		S <- createMatrixS(model)
	} else {
		S <- addVariablesAS(S, model, newLatent, newManifest)
	}
	if (!is.null(M)) {
		M <- addVariablesM(M, model, newLatent, newManifest)
	}
	model[['A']] <- A
	model[['S']] <- S
	model[['M']] <- M
	model[['F']] <- createMatrixF(model)
	return(model)
}


listArgumentRAM <- function(model, lst, remove) {
	if(remove == TRUE) {
		model <- removeEntriesRAM(model, lst)
	} else {
		model <- addEntriesRAM(model, lst)
	}
	return(model)
}

addEntriesRAM <- function(model, entries) {
	if (length(entries) == 0) {
		return(model)
	}
	filter <- sapply(entries, is, "MxPath")
	paths <- entries[filter]
	checkPaths(model, paths)
	if (length(paths) > 0) {
		for(i in 1:length(paths)) {
			model <- insertPathRAM(model, paths[[i]])
		}
	}
	filter <- sapply(entries, is, "MxData")
	data <- entries[filter]
	if (length(data) > 0) {
		if (length(data) > 1) {
			warning("Multiple data sources specified.  Only one will be chosen.")
		}
		data <- data[[1]]
		model$data <- data
	}
	return(model)
}

requireMeansVector <- function(data) {
	return(!is.null(data) && ((data@type == 'raw') ||
		((data@type == 'cov' || data@type == 'cor') &&
		 !(length(data@means) == 1 && is.na(data@means)))))
}

removeEntriesRAM <- function(model, entries) {
	if (length(entries) == 0) {
		return(model)
	}
	filter <- sapply(entries, is, "MxPath")
	paths <- entries[filter]
	if (length(paths) > 0) for(i in 1:length(paths)) {
		model <- removePathRAM(model, paths[[i]])
	}
	return(model)
}

getNotPathsOrData <- function(lst) {
	if (length(lst) == 0) {
		return(lst)
	}
	pathfilter <- sapply(lst, is, "MxPath")
	datafilter <- sapply(lst, is, "MxData")
	retval <- lst[!(pathfilter | datafilter)]
	return(retval)
}

checkPaths <- function(model, paths) {
	variables <- c(model@manifestVars, model@latentVars)
	fromNames <- unlist(sapply(paths, slot, 'from'))
	toNames <- unlist(sapply(paths, slot, 'to'))
	if(is.null(fromNames)) { fromNames <- character() }
	if(is.null(toNames)) { toNames <- character() }
	if (any(is.na(fromNames)) || any(is.na(toNames))) {
		stop("The \'from\' field or the \'to\' field contains an NA", call. = FALSE)
	}
	missingSource <- setdiff(fromNames, variables)
	missingSink   <- setdiff(toNames, variables)
	missingSource <- setdiff(missingSource, "one")
	if(length(missingSource) > 0) {
		stop(paste("The following are neither manifest nor latent variables:",
			omxQuotes(missingSource)), call. = FALSE)
	}
	if(length(missingSink) > 0) {
		stop(paste("The following are neither manifest nor latent variables:",
			omxQuotes(missingSink)), call. = FALSE)
	}
}

insertPathRAM <- function(model, path) {
	allfrom <- path@from
	allto <- path@to
	allarrows <- path@arrows
	excludeself <- path@excludeself
	if (length(allfrom) == 1 && (allfrom == "one")) {
		model <- insertMeansPathRAM(model, path)
		return(model)
	}
	default1 <- getOption('mxRAMDefaultSingle')
	default2 <- getOption('mxRAMDefaultDouble')
	A <- model[['A']]
	S <- model[['S']]
	F <- model[['F']]
	maxlength <- max(length(allfrom), length(allto))
	if (is.null(A)) { A <- createMatrixA(model) }
	if (is.null(S)) { S <- createMatrixS(model) }
	if (is.null(F)) { F <- createMatrixF(model) }
	for(i in 0:(maxlength - 1)) {
		from <- allfrom[i %% length(allfrom) + 1]
		to <- allto[i %% length(allto) + 1]
		arrows <- allarrows[i %% length(allarrows) + 1]
		if (excludeself && identical(from, to)) next
		if (arrows == 1) {
			A <- matrixSetPath(A, from, to, path, i, default1)
			S <- matrixClearPath(S, from, to)
			S <- matrixClearPath(S, to, from)
		} else if (arrows == 2) {
			S <- matrixSetPath(S, from, to, path, i, default2)
			S <- matrixSetPath(S, to, from, path, i, default2)
			A <- matrixClearPath(A, from, to)
			A <- matrixClearPath(A, to, from)
		} else {
			stop(paste("Unknown arrow type", arrows, 
					"with source", omxQuotes(from), 
				"and sink", omxQuotes(to)),
				call. = FALSE)
		}
	}
	model[['A']] <- A
	model[['S']] <- S
	model[['F']] <- F
	return(model)
}

insertMeansPathRAM <- function(model, path) {
	allto <- path@to
	arrows <- path@arrows
	if (any(arrows != 1)) {
		stop(paste('The means path to variable', omxQuotes(to),
			'does not contain a single-headed arrow.'), call. = FALSE)
	}
	M <- model$M
	if (is.null(M)) { 
		if(!is.null(model@objective) && is(model@objective,"MxRAMObjective") &&
			is.na(model@objective@M)) {
				model@objective@M <- "M"
		}
		M <- createMatrixM(model)
	}
	for(i in 0:(length(allto) - 1)) {
		to <- allto[i %% length(allto) + 1]
		M <- matrixSetPath(M, to, 1, path, i, 0)
	}
	model$M <- M
	return(model)
}

removePathRAM <- function(model, path) {
	allfrom <- path@from
	allto <- path@to
	allarrows <- path@arrows
	excludeself <- path@selfexlude
	A <- model[['A']]
	S <- model[['S']]
	F <- model[['F']]
	if (is.null(A)) { A <- createMatrixA(model) }
	if (is.null(S)) { S <- createMatrixS(model) }
	if (is.null(F)) { F <- createMatrixF(model) }
	maxlength <- max(length(allfrom), length(allto))
	for(i in 0:(maxlength - 1)) {
		from <- allfrom[i %% length(allfrom) + 1]
		to <- allto[i %% length(allto) + 1]
		arrows <- allarrows[i %% length(allarrows) + 1]
		if (excludeself && identical(from, to)) next
		if (arrows == 1) {
			A <- matrixClearPath(A, from, to)
		} else if (arrows == 2) {
			S <- matrixClearPath(S, from, to)
			S <- matrixClearPath(S, to, from)
		} else {
			stop(paste("Unknown arrow type", arrows, 
					"with source", omxQuotes(from), 
					"and sink", omxQuotes(to)),
					call. = FALSE)
		}
	}

	model[['A']] <- A
	model[['S']] <- S
	model[['F']] <- F
	return(model)
}

matrixSetPath <- function(mxMatrix, from, to, path, offset, default) {
	labels <- path@labels[offset %% length(path@labels) + 1]
	value <- path@values[offset %% length(path@values) + 1]
	free <- path@free[offset %% length(path@free) + 1]
	ubound <- path@ubound[offset %% length(path@ubound) + 1]
	lbound <- path@lbound[offset %% length(path@lbound) + 1]
	mxMatrix@values[to, from] <- value
	mxMatrix@labels[to, from] <- labels
	mxMatrix@ubound[to, from] <- ubound
	mxMatrix@lbound[to, from] <- lbound
	mxMatrix@free[to, from] <- free
	return(mxMatrix)
}

matrixClearPath <- function(mxMatrix, from, to) {
	mxMatrix@values[to, from] <- 0
	mxMatrix@labels[to, from] <- as.character(NA)
	mxMatrix@free[to, from] <- FALSE
	return(mxMatrix)
}

createMatrixM <- function(model) {
	variables <- c(model@manifestVars, model@latentVars)
	len <- length(variables)
	names <- list(NULL, variables)
	values <- matrix(0, 1, len)
	labels <- matrix(as.character(NA), 1, len)
	free <- matrix(c(rep.int(FALSE, length(model@manifestVars)),
		rep.int(FALSE, length(model@latentVars))), 1, len)
	retval <- mxMatrix("Full", values = values, free = free, labels = labels, name = "M")
	dimnames(retval) <- names
	return(retval)
}

createMatrixA <- function(model) {
	variables <- c(model@manifestVars, model@latentVars)
	len <- length(variables)
	names <- list(variables, variables)
	values <- matrix(0, len, len)
	free <- matrix(FALSE, len, len)
	labels <- matrix(as.character(NA), len, len)
	retval <- mxMatrix("Full", values = values, free = free, labels = labels, name = "A")
	dimnames(retval) <- names
	return(retval)
}

createMatrixS <- function(model) {
	variables <- c(model@manifestVars, model@latentVars)
	len <- length(variables)
	names <- list(variables, variables)
	values <- matrix(0, len, len)
	free <- matrix(FALSE, len, len)
	labels <- matrix(as.character(NA), len, len)
	retval <- mxMatrix("Symm", values = values, free = free, labels = labels, name = "S")
	dimnames(retval) <- names
	return(retval)
}

createMatrixF <- function(model) {
	variables <- c(model@manifestVars, model@latentVars)
	len <- length(variables)
	names <- list(model@manifestVars, variables)
	values <- diag(nrow = length(model@manifestVars), ncol = len)
	free <- matrix(FALSE, length(model@manifestVars), len)
	labels <- matrix(as.character(NA), length(model@manifestVars), len)
	retval <- mxMatrix("Full", values = values, free = free, labels = labels, name = "F")
	dimnames(retval) <- names
	return(retval)
}

addVariablesAS <- function(oldmatrix, model, newLatent, newManifest) {
	newLatent <- length(newLatent)
	newManifest <- length(newManifest)
	oldmatrix@values <- addVariablesMatrix(oldmatrix@values, 0, 
		model, newLatent, newManifest)
	oldmatrix@free <- addVariablesMatrix(oldmatrix@free, FALSE, 
		model, newLatent, newManifest)
	oldmatrix@labels <- addVariablesMatrix(oldmatrix@labels, as.character(NA), 
		model, newLatent, newManifest)
	oldmatrix@lbound <- addVariablesMatrix(oldmatrix@lbound, as.numeric(NA), 
		model, newLatent, newManifest)
	oldmatrix@ubound <- addVariablesMatrix(oldmatrix@ubound, as.numeric(NA), 
		model, newLatent, newManifest)		
	variables <- c(model@manifestVars, model@latentVars)
	dimnames(oldmatrix) <- list(variables, variables)
	return(oldmatrix)
}

addVariablesM <- function(oldmatrix, model, newLatent, newManifest) {
	oldmatrix@values <- addVariablesMatrixM(oldmatrix@values, 0, 0, model, newLatent, newManifest)
	oldmatrix@free   <- addVariablesMatrixM(oldmatrix@free, FALSE, TRUE, model, newLatent, newManifest)
	oldmatrix@labels <- addVariablesMatrixM(oldmatrix@labels, as.character(NA), as.character(NA),
		model, newLatent, newManifest) 
	oldmatrix@lbound <- addVariablesMatrixM(oldmatrix@lbound, as.numeric(NA), as.numeric(NA), 
		model, newLatent, newManifest)
	oldmatrix@ubound <- addVariablesMatrixM(oldmatrix@ubound, as.numeric(NA), as.numeric(NA), 
		model, newLatent, newManifest)
	dimnames(oldmatrix) <- list(NULL, c(model@manifestVars, model@latentVars))
	return(oldmatrix)
}

removeVariablesAS <- function(oldmatrix, variables) {
	if (length(variables) > 0) {
		for (i in 1:length(variables)) {
			index <- match(variables[[i]], dimnames(oldmatrix)[[1]])
			oldmatrix@values <- oldmatrix@values[-index, -index]
			oldmatrix@free <- oldmatrix@free[-index, -index]
			oldmatrix@labels <- oldmatrix@labels[-index, -index]
			oldmatrix@lbound <- oldmatrix@lbound[-index, -index]
			oldmatrix@ubound <- oldmatrix@ubound[-index, -index]
		}
	}
	return(oldmatrix)
}


addVariablesMatrix <- function(oldmatrix, value, model, newLatent, newManifest) {
	currentManifest <- length(model@manifestVars) - newManifest
	currentLatent <- length(model@latentVars) - newLatent
	newSize <- length(model@manifestVars) + length(model@latentVars)
	if (currentManifest > 0) {
		manifestXmanifest <- oldmatrix[1 : currentManifest, 1 : currentManifest]
	} else {
		manifestXmanifest <- matrix(value, 0, 0)
	}
	if (currentLatent > 0) {
		latentStart <- currentManifest + 1
		latentEnd <- currentManifest + currentLatent
		manifestXlatent <- oldmatrix[1 : currentManifest, latentStart : latentEnd]
		latentXmanifest <- oldmatrix[latentStart : latentEnd, 1 : currentManifest]
		latentXlatent <- oldmatrix[latentStart : latentEnd, latentStart : latentEnd]
	} else {
		manifestXlatent <- matrix(value, 0, 0)
		latentXmanifest <- matrix(value, 0, 0)
		latentXlatent <- matrix(value, 0, 0)
	}
	newtop <- cbind(manifestXmanifest, matrix(value, currentManifest, newManifest),
					manifestXlatent, matrix(value, currentManifest, newLatent)) 
	newtop <- rbind(newtop, matrix(value, newManifest, newSize))
	newbottom <- cbind(latentXmanifest, matrix(value, currentLatent, newManifest),
					latentXlatent, matrix(value, currentLatent, newLatent)) 
	newbottom <- rbind(newbottom, matrix(value, newLatent, newSize))
	newmatrix <- rbind(newtop, newbottom)
	return(newmatrix)
}

addVariablesMatrixM <- function(oldmatrix, newLatentValue, newManifestValue, model, newLatent, newManifest) {
	newManifest <- length(newManifest)
	newLatent <- length(newLatent)
	currentManifest <- length(model@manifestVars) - newManifest
	currentLatent <- length(model@latentVars) - newLatent
	values <- c(oldmatrix[1, 1:currentManifest], 
		rep.int(newManifestValue, newManifest),
		oldmatrix[1, (currentManifest + 1) : (currentLatent + currentManifest)],
		rep.int(newLatentValue, newLatent))
	newmatrix <- matrix(values, 1, length(model@manifestVars) + length(model@latentVars))
	return(newmatrix)
}

replaceMethodRAM <- function(model, index, value) {
	pair <- omxReverseIdentifier(model, index)
	namespace <- pair[[1]]
	name <- pair[[2]]
	if (namespace == model@name && name == "data") {
		model@data <- value
		if (requireMeansVector(value)) {
			model@objective@M <- "M"
		}
	} else {
		model <- omxReplaceMethod(model, index, value)
	}
	return(model)
}
