#
#   Copyright 2007-2013 The OpenMx Project
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

imxModelTypes[['RAM']] <- "MxRAMModel"

imxVariableTypes <- c(imxVariableTypes, "exogenous", "endogenous")

# Define generic functions

setMethod("imxTypeName", "MxRAMModel", 
	function(model) { "RAM" }
)

setMethod("imxInitModel", "MxRAMModel", 
	function(model) {
		if (is.null(model[['A']])) {
			model[['A']] <- createMatrixA(model)
		}
		if (is.null(model[['S']])) {
			model[['S']] <- createMatrixS(model)
		}
		if (is.null(model[['expectation']])) {
			model[['expectation']] <- mxExpectationRAM('A', 'S', 'F')
		}
		if (is.null(model[['fitfunction']])) {
			model[['fitfunction']] <- mxFitFunctionML()
		}
		model[['F']] <- createMatrixF(model)
		return(model)
	}
)

setMethod("imxModelBuilder", "MxRAMModel", 
	function(model, lst, name, 
		manifestVars, latentVars, remove, independent) {
        model <- nameArgument(model, name)		
		model <- variablesArgumentRAM(model, manifestVars, latentVars, remove)
		model <- listArgumentRAM(model, lst, remove)
		notPathOrData <- getNotPathsOrData(lst)
		callNextMethod(model, notPathOrData, NA, character(), 
			character(), remove, independent)
	}
)

setMethod("imxVerifyModel", "MxRAMModel",
	function(model) {
		if ((length(model$A) == 0) ||
			(length(model$S) == 0) ||
			(length(model$F) == 0)) {
				msg <- paste("The RAM model", omxQuotes(model@name),
                "does not contain any paths.",
				" You can add paths to your model like this:",
				" mxPath(from = 'x1', to = 'y1')")
				stop(msg, call. = FALSE)
		}
		expectation <- model$expectation
		if (!is.null(expectation) && is(expectation, "MxExpectationRAM")) {
			if (!is.null(model@data) && model@data@type == "raw" &&
   	    	is.null(model$M)) {
				msg <- paste("The RAM model", omxQuotes(model@name),
        	       "contains raw data but has not specified any means paths.",
                   "Add something like mxPath(from = 'one', to = manifests) to your model."
                   )
				stop(msg, call. = FALSE)
			}
			if (!is.null(model@data) && !single.na(model@data@means) &&
				is.null(model$M)) {
				msg <- paste("The RAM model", omxQuotes(model@name),
					"contains an observed means vector",
					"but has not specified any means paths.")
				stop(msg, call. = FALSE)        	
			}
		}
		if (length(model@submodels) > 0) {
			return(all(sapply(model@submodels, imxVerifyModel)))
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
		latentVars <- varsToCharacter(latentVars, "latent")
		manifestVars <- varsToCharacter(manifestVars, "manifest")
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

	modelLatent   <- unlist(model@latentVars, use.names = FALSE)
	modelManifest <- unlist(model@manifestVars, use.names = FALSE)

	model <- addVariablesHelper(model, "latentVars", latent)
	model <- addVariablesHelper(model, "manifestVars", manifest)

	latent <- unlist(latent, use.names = FALSE)
	manifest <- unlist(manifest, use.names = FALSE)

	newLatent   <- setdiff(latent, modelLatent)
	newManifest <- setdiff(manifest, modelManifest)

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
	if (length(paths) > 0) {
		model <- insertAllPathsRAM(model, paths)
	}
	filter <- sapply(entries, is, "MxData")
	data <- entries[filter]
	if (length(data) > 0) {
		if (length(data) > 1) {
			warning("Multiple data sources specified.  Only one will be chosen.")
		}
		data <- data[[1]]
		model@data <- data
		model[['F']] <- createMatrixF(model)
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
	if (length(paths) > 0) {
		model <- removeAllPathsRAM(model, paths)
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

expectationIsMissingMeans <- function(model) {
	expectation <- model@expectation
	return(!is.null(expectation) &&
		is(expectation, "MxExpectationRAM") &&
		is.na(expectation@M))
}

insertAllPathsRAM <- function(model, paths) {
	A <- model[['A']]
	S <- model[['S']]
	M <- model[['M']]
	if (is.null(A)) { A <- createMatrixA(model) }
	if (is.null(S)) { S <- createMatrixS(model) }
	
	legalVars <- c(model@latentVars, model@manifestVars, "one")
	
	for(i in 1:length(paths)) {
		path <- paths[[i]]
	
		missingvalues <- is.na(path@values)
		path@values[missingvalues] <- 0
		
		if (single.na(path@to)) {
			path@to <- path@from
			paths[[i]] <- path
		}
		
		allFromTo <- unique(c(path@from, path@to))
                varExist <- allFromTo %in% legalVars 
		if(!all(varExist)) {
			missingVars <- allFromTo[!varExist]
			stop(paste("Nice try, you need to add", 
				omxQuotes(missingVars), 
				"to either manifestVars or latentVars before you",
				"can use them in a path."), call. = FALSE)
		}
		
		if (length(path@from) == 1 && (path@from == "one")) {
			if (is.null(M)) {
				M <- createMatrixM(model) 
				if(expectationIsMissingMeans(model)) {
					model@expectation@M <- "M"
				}
			}
			M <- insertMeansPathRAM(path, M)
		} else {
			expanded <- expandPathConnect(path@from, path@to, path@connect)
			path@from <- expanded$from
			path@to   <- expanded$to
			retval <- insertPathRAM(path, A, S)
			A <- retval[[1]]
			S <- retval[[2]]	
		}
	}
	checkPaths(model, paths)
	model[['A']] <- A
	model[['S']] <- S
	if (!is.null(M)) {
		model[['M']] <- M
	}
	
	return(model)
}

removeAllPathsRAM <- function(model, paths) {
	A <- model[['A']]
	S <- model[['S']]
	M <- model[['M']]
	if (is.null(A)) { A <- createMatrixA(model) }
	if (is.null(S)) { S <- createMatrixS(model) }
	for(i in 1:length(paths)) {

		path <- paths[[i]]

		if (single.na(path@to)) {
			path@to <- path@from
			paths[[i]] <- path
		}
		
		if (length(path@from) == 1 && (path@from == "one")) {		
			M <- removeMeansPathRAM(path, M)
		} else {
			expanded <- expandPathConnect(path@from, path@to, path@connect)
			path@from <- expanded$from
			path@to   <- expanded$to
			retval <- removePathRAM(path, A, S)
			A <- retval[[1]]
			S <- retval[[2]]
		}
	}
	checkPaths(model, paths)
	model[['A']] <- A
	model[['S']] <- S
	if (!is.null(M)) {
		model[['M']] <- M
	}
	return(model)
}


insertPathRAM <- function(path, A, S) {
	allfrom <- path@from
	allto <- path@to
	allarrows <- path@arrows
	allfree <- path@free
	allvalues <- path@values
	alllabels <- path@labels
	alllbound <- path@lbound
	allubound <- path@ubound
	maxlength <- max(length(allfrom), length(allto))
	A_free <- A@free
	A_values <- A@values
	A_labels <- A@labels
	A_lbound <- A@lbound
	A_ubound <- A@ubound
	S_free   <- S@free
	S_values <- S@values
	S_labels <- S@labels
	S_lbound <- S@lbound
	S_ubound <- S@ubound
	for(i in 0:(maxlength - 1)) {
		from <- allfrom[[i %% length(allfrom) + 1]]
		to <- allto[[i %% length(allto) + 1]]
		arrows <- allarrows[[i %% length(allarrows) + 1]]
		nextvalue <- allvalues[[i %% length(allvalues) + 1]]
		nextfree <- allfree[[i %% length(allfree) + 1]]
		nextlabel <- alllabels[[i %% length(alllabels) + 1]]
		nextubound <- allubound[[i %% length(allubound) + 1]]
		nextlbound <- alllbound[[i %% length(alllbound) + 1]]		
		if (arrows == 1) {
			A_free[to, from] <- nextfree
			A_values[to, from] <- nextvalue
			A_labels[to, from] <- nextlabel
			A_ubound[to, from] <- nextubound
			A_lbound[to, from] <- nextlbound
			S_values[to, from] <- 0
			S_labels[to, from] <- as.character(NA)
			S_free[to, from] <- FALSE
			S_values[from, to] <- 0
			S_labels[from, to] <- as.character(NA)
			S_free[from, to] <- FALSE			
		} else if (arrows == 2) {
			S_free[to, from] <- nextfree
			S_values[to, from] <- nextvalue
			S_labels[to, from] <- nextlabel
			S_ubound[to, from] <- nextubound
			S_lbound[to, from] <- nextlbound
			S_free[from, to] <- nextfree
			S_values[from, to] <- nextvalue
			S_labels[from, to] <- nextlabel
			S_ubound[from, to] <- nextubound
			S_lbound[from, to] <- nextlbound
			A_values[to, from] <- 0
			A_labels[to, from] <- as.character(NA)
			A_free[to, from] <- FALSE
			A_values[to, from] <- 0
			A_labels[to, from] <- as.character(NA)
			A_free[to, from] <- FALSE
		} else {
			stop(paste("Unknown arrow type", arrows, 
					"with source", omxQuotes(from), 
				"and sink", omxQuotes(to)),
				call. = FALSE)
		}
	}
	A@free <- A_free
	A@values <-	A_values 
	A@labels <-	A_labels 
	A@lbound <-	A_lbound 
	A@ubound <-	A_ubound 
	S@free <- S_free   
	S@values <-	S_values 
	S@labels <-	S_labels 
	S@lbound <- S_lbound 
	S@ubound <- S_ubound 
	return(list(A, S))
}

removeMeansPathRAM <- function(path, M) {
	if(is.null(M)) {
		return(NULL)
	}
	allto <- path@to
	for(i in 0:(length(allto) - 1)) {
		to <- allto[[i %% length(allto) + 1]]
		M@free[1, to] <- FALSE
		M@values[1, to] <- 0
		M@labels[1, to] <- as.character(NA)
	}
	return(M)
}

insertMeansPathRAM <- function(path, M) {
	allto     <- path@to
	allarrows <- path@arrows
	allfree   <- path@free
	allvalues <- path@values
	alllabels <- path@labels
	alllbound <- path@lbound
	allubound <- path@ubound	
	if (any(allarrows != 1)) {
		stop(paste('The means path must be a single-headed arrow\n',
		'path from "one" to variable(s)', omxQuotes(allto)), call. = FALSE)
	}
	for(i in 0:(length(allto) - 1)) {
		to <- allto[[i %% length(allto) + 1]]
		nextvalue  <- allvalues[[i %% length(allvalues) + 1]]
		nextfree   <- allfree[[i %% length(allfree) + 1]]
		nextlabel  <- alllabels[[i %% length(alllabels) + 1]]
		nextubound <- allubound[[i %% length(allubound) + 1]]
		nextlbound <- alllbound[[i %% length(alllbound) + 1]]
		M@free[1, to] <- nextfree
		M@values[1, to] <- nextvalue
		M@labels[1, to] <- nextlabel
		M@ubound[1, to] <- nextubound
		M@lbound[1, to] <- nextlbound
	}
	return(M)
}

removePathRAM <- function(path, A, S) {
	allfrom <- path@from
	allto <- path@to
	allarrows <- path@arrows
	maxlength <- max(length(allfrom), length(allto))
	A_free <- A@free
	A_values <- A@values
	A_labels <- A@labels
	S_free   <- S@free
	S_values <- S@values
	S_labels <- S@labels
	for(i in 0:(maxlength - 1)) {
		from <- allfrom[i %% length(allfrom) + 1]
		to <- allto[i %% length(allto) + 1]
		arrows <- allarrows[i %% length(allarrows) + 1]
		if (arrows == 1) {
			A_values[to, from] <- 0
			A_labels[to, from] <- as.character(NA)
			A_free[to, from] <- FALSE		
		} else if (arrows == 2) {
			S_values[to, from] <- 0
			S_labels[to, from] <- as.character(NA)
			S_free[to, from] <- FALSE		
			S_values[from, to] <- 0
			S_labels[from, to] <- as.character(NA)
			S_free[from, to] <- FALSE					
		} else {
			stop(paste("Unknown arrow type", arrows, 
					"with source", omxQuotes(from), 
					"and sink", omxQuotes(to)),
					call. = FALSE)
		}
	}
	A@free <- A_free
	A@values <-	A_values 
	A@labels <-	A_labels 
	S@free <- S_free   
	S@values <-	S_values 
	S@labels <-	S_labels 	
	return(list(A, S))
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
	values <- diag(nrow = length(model@manifestVars), ncol = len)
	names <- list(model@manifestVars, variables)
	if (!is.null(model@data) && (model@data@type != 'raw')) {
		manifestNames <- rownames(model@data@observed)
		extraData <- setdiff(manifestNames, model@manifestVars)
		extraVars <- setdiff(model@manifestVars, manifestNames)
		if (length(extraData) > 0) {
			msg <- paste("The observed data contains the variables:",
				omxQuotes(extraData), "that have not been declared in the",
				"manifest variables.")
			stop(msg, call. = FALSE)
		}
		if (length(extraVars) > 0) {
			msg <- paste("The manifest variables include",
				omxQuotes(extraVars), "that have not been found in the",
				"observed data.")
			stop(msg, call. = FALSE)
		}
		dimnames(values) <- names
		values <- values[manifestNames,]
		names <- list(manifestNames, variables)
	}
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
	pair <- imxReverseIdentifier(model, index)
	namespace <- pair[[1]]
	name <- pair[[2]]
	if (namespace == model@name && name == "data") {
		model@data <- value
		if (requireMeansVector(value)) {
			model@expectation@M <- "M"
		}
	} else {
		model <- imxReplaceMethod(model, index, value)
	}
	return(model)
}

