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
		if (is.null(model[['objective']])) {
			model[['objective']] <- mxRAMObjective()
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
	if (!is.null(model[['M']])) {
		model[['M']] <- addVariablesM(model[['M']], 
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
	filter <- sapply(entries, omxIsPath)
	paths <- entries[filter]
	checkPaths(model, paths)
	if (length(paths) > 0) for(i in 1:length(paths)) {
		model <- insertPathRAM(model, paths[[i]])
	}
	filter <- sapply(entries, function(x) { is(x, "MxData") })
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

useMeansVector <- function(data) {
	return(!is.null(data) && ((data@type == 'raw') ||
		((data@type == 'cov' || data@type == 'cor') &&
		 !(length(data@means) == 1 && is.na(data@means)))))
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

getNotPathsOrData <- function(lst) {
	filter <- sapply(lst, function(x) {
		omxIsPath(x) || is(x, "MxData")})
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
	from <- path$from
	to <- path$to
	arrows <- path$arrows
	if (from == "one") {
		model <- insertMeansPathRAM(model, path)
		return(model)
	}
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

insertMeansPathRAM <- function(model, path) {
	to <- path$to
	arrows <- path$arrows
	if (arrows != 1) {
		stop(paste('The means path to variable', omxQuotes(to),
			'does not contain a single-headed arrow.'), call. = FALSE)
	}
	if (is.null(model[['M']])) { model[['M']] <- createMatrixM(model) }
	model[['M']] <- matrixSetPath(model[['M']], to, 1, path, 0)
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
	labels <- path$labels
	value <- path$values
	free <- path$free
	ubound <- path$ubound
	lbound <- path$lbound
	if (is.null(value)) {
		mxMatrix@values[to, from] <- as.numeric(default)
	} else {
		mxMatrix@values[to, from] <- as.numeric(value)
	}
	if (is.null(labels)) {
		mxMatrix@labels[to, from] <- as.character(NA)
	} else {
		mxMatrix@labels[to, from] <- as.character(labels)
	}
	if (!is.null(ubound)) {
		mxMatrix@ubound[to, from] <- as.numeric(ubound)
	}
	if (!is.null(lbound)) {
		mxMatrix@lbound[to, from] <- as.numeric(lbound)
	}
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
	free <- matrix(c(rep.int(TRUE, length(model@manifestVars)),
		rep.int(FALSE, length(model@latentVars))), 1, len)
	retval <- mxMatrix("Full", values, free, labels, name = "M")
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
	retval <- mxMatrix("Full", values, free, labels, name = "A")
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
	retval <- mxMatrix("Symm", values, free, labels, name = "S")
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
	retval <- mxMatrix("Full", values, free, labels, name = "F")
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
	values <- c(oldmatrix[1:currentManifest, 1], 
		rep.int(newManifestValue, newManifest),
		oldmatrix[(currentManifest + 1) : (currentLatent + currentManifest), 1],
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
		if (!useMeansVector(value)) {
			model[['M']] <- NULL
			if(!is.null(model@objective) && is(model@objective,"MxRAMObjective") &&
				!is.na(model@objective@M)) {
					model@objective@M <- as.character(NA)
			}
		} else if (is.null(model[['M']])) {
			model[['M']] <- createMatrixM(model)
			if(!is.null(model@objective) && is(model@objective,"MxRAMObjective") &&
				is.na(model@objective@M)) {
					model@objective@M <- "M"
			}
		}
	} else {
		model <- omxReplaceMethod(model, index, value)
	}
	return(model)
}