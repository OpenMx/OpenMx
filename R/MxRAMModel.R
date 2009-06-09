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
		model[['F']] <- createMatrixF(model)
		return(model)
	}
)

setMethod("omxModelBuilder", "MxRAMModel", 
	function(model, lst, name, 
		manifestVars, latentVars, remove, independent) {
		model <- variablesArgumentRAM(model, manifestVars, latentVars, remove)
		model <- listArgumentRAM(model, lst, remove)
		notpaths <- getNotPathsRAM(lst)
		callNextMethod(model, notpaths, name, character(), 
			character(), remove, independent)
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
	model[['A']] <- removeVariablesAS(model[['A']], latent)
	model[['A']] <- removeVariablesAS(model[['A']], manifest)
	model[['S']] <- removeVariablesAS(model[['S']], latent)
	model[['S']] <- removeVariablesAS(model[['S']], manifest)
	model[['F']] <- createMatrixF(model)
	return(model)
}

addVariablesRAM <- function(model, latent, manifest) {
	newLatent   <- setdiff(latent, model@latentVars)
	newManifest <- setdiff(manifest, model@manifestVars)
	model@latentVars <- c(model@latentVars, newLatent)
	model@manifestVars <- c(model@manifestVars, newManifest)
	if (is.null(model[['A']])) {
		model[['A']] <- createMatrixA(model)
	} else {
		model[['A']] <- addVariablesAS(model[['A']], 
			model, newLatent, newManifest)
	}
	if (is.null(model[['S']])) {
		model[['S']] <- createMatrixS(model)
	} else {
		model[['S']] <- addVariablesAS(model[['S']], 
			model, newLatent, newManifest)
	}
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
	filter   <- sapply(entries, omxIsPath)
	paths    <- entries[filter]
	checkPaths(model, paths)
	if (length(paths) > 0) for(i in 1:length(paths)) {
		model <- insertPathRAM(model, paths[[i]])
	}
	return(model)
}

removeEntriesRAM <- function(model, entries) {
	if (length(entries) == 0) {
		return(model)
	}
	filter   <- sapply(entries, omxIsPath)
	paths    <- entries[filter]
	notpaths <- entries[!filter]
	if (length(paths) > 0) for(i in 1:length(paths)) {
		model <- removePathRAM(model, paths[[i]])
	}
	return(model)
}

getNotPathsRAM <- function(lst) {
	filter <- sapply(lst, omxIsPath)
	notpaths <- lst[!filter]
	return(notpaths)
}

pathSelectFrom <- function(lst) {
  retval <- lapply(lst, function(x) { x$from } )
  return(retval)
}	

pathSelectTo <- function(lst) {
  retval <- lapply(lst, function(x) { x$to } )
  return(retval)
}

checkPaths <- function(model, paths) {
	variables <- c(model@manifestVars, model@latentVars)
	if (any(is.na(pathSelectFrom(paths))) || any(is.na(pathSelectTo(paths)))) {
		stop("The \'from\' field or the \'to\' field contains an NA", call. = FALSE)
	}
	missingSource <- setdiff(pathSelectFrom(paths), variables)
	missingSink   <- setdiff(pathSelectTo(paths), variables)
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
	from <- path$from
	to <- path$to
	arrows <- path$arrows
	default1 <- getOption('mxRAMDefaultSingle')
	default2 <- getOption('mxRAMDefaultDouble')
	if (is.null(model[['A']])) { model[['A']] <- createMatrixA(model) }
	if (is.null(model[['S']])) { model[['S']] <- createMatrixS(model) }
	if (is.null(model[['F']])) { model[['F']] <- createMatrixF(model) }
	if (arrows == 1) {
		model[['A']] <- matrixSetPath(model[['A']], from, to, path, default1)
		model[['S']] <- matrixClearPath(model[['S']], from, to)
		model[['S']] <- matrixClearPath(model[['S']], to, from)
	} else if (arrows == 2) {
		model[['S']] <- matrixSetPath(model[['S']], from, to, path, default2)
		model[['S']] <- matrixSetPath(model[['S']], to, from, path, default2)
		model[['A']] <- matrixClearPath(model[['A']], from, to)
		model[['A']] <- matrixClearPath(model[['A']], to, from)
	} else {
		stop(paste("Unknown arrow type", arrows, 
				"with source", omxQuotes(from), 
				"and sink", omxQuotes(to)),
				call. = FALSE)
	}
	return(model)
}

removePathRAM <- function(model, path) {
	from <- path$from
	to <- path$to
	arrows <- path$arrows
	if (is.null(model[['A']])) { model[['A']] <- createMatrixA(model) }
	if (is.null(model[['S']])) { model[['S']] <- createMatrixS(model) }
	if (is.null(model[['F']])) { model[['F']] <- createMatrixF(model) }
	if (arrows == 1) {
		model[['A']] <- matrixClearPath(model[['A']], from, to)
	} else if (arrows == 2) {
		model[['S']] <- matrixClearPath(model[['S']], from, to)
		model[['S']] <- matrixClearPath(model[['S']], to, from)
	} else {
		stop(paste("Unknown arrow type", arrows, 
				"with source", omxQuotes(from), 
				"and sink", omxQuotes(to)),
				call. = FALSE)
	}
	return(model)
}

matrixSetPath <- function(mxMatrix, from, to, path, default) {
	name <- path$name
	value <- path$start
	free <- path$free
	if (is.null(value)) {
		mxMatrix@values[to, from] <- default		
	} else {
		mxMatrix@values[to, from] <- value
	}
	if (is.null(name)) {
		mxMatrix@labels[to, from] <- NA
	} else {
		mxMatrix@labels[to, from] <- name
	}
	mxMatrix@free[to, from] <- free
	return(mxMatrix)
}

matrixClearPath <- function(mxMatrix, from, to) {
	mxMatrix@values[to, from] <- 0
	mxMatrix@labels[to, from] <- NA
	mxMatrix@free[to, from] <- FALSE
	return(mxMatrix)
}

createMatrixA <- function(model) {
	variables <- c(model@manifestVars, model@latentVars)
	len <- length(variables)
	names <- list(variables, variables)
	values <- matrix(0, len, len, dimnames = names)
	free <- matrix(FALSE, len, len, dimnames = names)
	labels <- matrix(as.character(NA), len, len, dimnames = names)
	retval <- mxMatrix("Full", values, free, labels, name = "A")
	return(retval)
}

createMatrixS <- function(model) {
	variables <- c(model@manifestVars, model@latentVars)
	len <- length(variables)
	names <- list(variables, variables)
	values <- matrix(0, len, len, dimnames = names)	
	free <- matrix(FALSE, len, len, dimnames = names)
	labels <- matrix(as.character(NA), len, len, dimnames = names)
	retval <- mxMatrix("Symm", values, free, labels, name = "S")
	return(retval)
}

createMatrixF <- function(model) {
	variables <- c(model@manifestVars, model@latentVars)
	len <- length(variables)
	names <- list(model@manifestVars, variables)
	values <- diag(nrow = length(model@manifestVars), ncol = len)
	free <- matrix(FALSE, length(model@manifestVars), len, dimnames = names)
	labels <- matrix(as.character(NA), length(model@manifestVars), len, dimnames = names)
	dimnames(values) <- names
	retval <- mxMatrix("Full", values, free, labels, name = "F")
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
	return(oldmatrix)
}

removeVariablesAS <- function(oldmatrix, variables) {	
	if (length(variables) > 0) {
		for (i in 1:length(variables)) {
			index <- match(variables[[i]], dimnames(oldmatrix@values)[[1]])
			oldmatrix@values <- oldmatrix@values[-index, -index]
			oldmatrix@free <- oldmatrix@free[-index, -index]
			oldmatrix@labels <- oldmatrix@labels[-index, -index]
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
	variables <- c(model@manifestVars, model@latentVars)
	dimnames(newmatrix) <- list(variables, variables)
	return(newmatrix)
}

