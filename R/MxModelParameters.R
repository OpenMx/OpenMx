#
#   Copyright 2007-2021 by the individuals mentioned in the source code history
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


omxLocateParameters <- function(model, labels = NULL, indep = FALSE,
                                free = c(TRUE, FALSE, NA)) {
	if (identical(free, c(TRUE, FALSE, NA))) {
		free <- TRUE
	}
	retval <- locateParametersInternal(model, labels, indep, free)
	if (nrow(retval) > 0) {
		retval <- retval[order(retval$label),]
	}
	rownames(retval) <- NULL
	return(retval)
}

locateParametersInternal <- function(model, labels, indep, free) {
	if (!is.null(labels) && !single.na(labels) && !is.character(labels)) {
		stop("'labels' argument must be NULL or a character vector")
	}
	retval <- lapply(model@matrices, locateParametersHelper, model@name, labels, free)
	retval <- do.call(rbind, retval)
	if(indep) {
		submodels <- model@submodels
	} else {
		submodels <- imxDependentModels(model)
	}
	if (length(submodels) > 0) {
		subparams <- lapply(submodels, locateParametersInternal, labels, indep, free)
		subparams <- do.call(rbind, subparams)
		retval <- rbind(retval, subparams)
	}
	return(retval)
}

locateParametersHelper <- function(matrix, modelname, target, free) {
	retval <- data.frame(label = character(0), model = character(0),
		matrix = character(0), row = numeric(0),
		col = numeric(0), value = numeric(0),
		lbound = numeric(0), ubound = numeric(0),
		stringsAsFactors = FALSE)
	freeMask <- matrix@free
  if (is.na(free)) freeMask[,] <- TRUE
  else freeMask <- freeMask == free
	count <- sum(freeMask)
	if (count == 0) {
		return(retval)
	} else {
		labels <- matrix@labels[freeMask]
		rows <- row(freeMask)[freeMask]
		cols <- col(freeMask)[freeMask]
		values <- matrix@values[freeMask]
		lbound <- matrix@lbound[freeMask]
		ubound <- matrix@ubound[freeMask]
		for(i in 1:count) {
			pname <- labels[[i]]
			if (is.null(target) || (is.na(pname) && any(is.na(target))) ||
				pname %in% target) {
				nextentry <- nrow(retval) + 1
				retval[nextentry,'label'] <- pname
				retval[nextentry,'model'] <- modelname
				retval[nextentry,'matrix'] <- matrix@name
				retval[nextentry,'row'] <- rows[[i]]
				retval[nextentry,'col'] <- cols[[i]]
				retval[nextentry,'value'] <- values[[i]]
				retval[nextentry,'lbound'] <- lbound[[i]]
				retval[nextentry,'ubound'] <- ubound[[i]]
			}
		}
	}
	return(retval)
}

omxGetParameters <- function(model, indep = FALSE, free = c(TRUE, FALSE, NA),
			     fetch = c('values', 'free', 'lbound', 'ubound', 'all'),
			     labels=c()) {
	if (identical(free, c(TRUE, FALSE, NA))) {
		free <- TRUE
	}
	if (identical(fetch, c('values', 'free', 'lbound', 'ubound', 'all'))) {
		fetch <- 'values'
	}
	if (!is.logical(free) || length(free) != 1) {
		stop("argument 'free' must be a 'TRUE', 'FALSE', or NA")
	}
	if (!is.character(fetch) || length(fetch) != 1 ||
		!(fetch %in% c('values', 'free', 'lbound', 'ubound', 'all'))) {
		stop("argument 'fetch' must be one of c('values', 'free', 'lbound', 'ubound', 'all')")
	}
	if (fetch == 'all') {
		values <- omxGetParameters(model, indep, free, 'values', labels)
		lbound <- omxGetParameters(model, indep, free, 'lbound', labels)
		ubound <- omxGetParameters(model, indep, free, 'ubound', labels)
		if (!is.na(free) && free) {
			free <- rep.int(TRUE, length(values))
		} else if (!is.na(free) && !free) {
			free <- rep.int(FALSE, length(values))
		} else {
			free <- omxGetParameters(model, indep, free, 'free')
		}
		return(data.frame(values, free, lbound, ubound))
	}
	parameters <- lapply(names(model@matrices), getParametersHelper, model, free, fetch, labels)
	plen <- lapply(parameters, length)
	parameters[plen == 0] <- NULL
	names(parameters) <- NULL
	parameters <- unlist(parameters)
	if(indep) {
		submodels <- model@submodels
	} else {
		submodels <- imxDependentModels(model)
	}
	if (length(submodels) > 0) {
		subparams <- lapply(submodels, omxGetParameters, indep, free, fetch, labels)
		plen <- lapply(subparams, length)
		subparams[plen == 0] <- NULL
		names(subparams) <- NULL
		subparams <- unlist(subparams)
		parameters <- c(parameters, subparams)
	}
	if (length(parameters) == 0) return(parameters)
	notDup1 <- !duplicated(names(parameters), incomparables = NA)
	notDup2 <- !duplicated(names(parameters[!is.na(parameters)]), incomparables = NA)
	if (sum(notDup1) == sum(notDup2)) {
		# Filter out NAs when some instance of the parameter does have a value.
		parameters <- parameters[!is.na(parameters)]
		notDup1 <- notDup2
	}
	parameters <- parameters[notDup1]
	return(parameters)
}

setParametersCheckVector <- function(values, test, argname, typename) {
	if (is.null(values)) return()
	if (!test(values)) {
		stop(paste(omxQuotes(argname),
			"argument must either be NA or a",
			typename, "vector"), call. = FALSE)
	}
}


omxSetParameters <- function(model, labels=names(coef(model)), free = NULL, values = NULL,
	newlabels = NULL, lbound = NULL, ubound = NULL, indep = FALSE,
	strict = TRUE, name = NULL) {
	if (missing(labels) || !is.character(labels) || length(labels) == 0) {
		stop("'labels' argument must be a character vector")
	}
	if (any(is.na(labels))) {
		stop("'labels' argument must not contain NA values")
	}
	if (any(duplicated(labels))) {
		stop("'labels' argument must not contain duplicate values")
	}
	if (!is.null(name) && length(name) != 1 && !is.character(name)) {
		stop("'name' argument must be a character string")
	}
	if (missing(free) && missing(values) && missing(newlabels) &&
		    missing(lbound) && missing(ubound)) {
		what <- c('free','values','newlabels','lbound','ubound')
		warning(paste("What do you want to change?",
			"Pick some of", omxQuotes(what)))
	}
	if (strict) {
		pnames <- names(omxGetParameters(model, indep, NA))
		missing <- setdiff(labels, pnames)
		if (length(missing) > 0) {
			msg <- paste("The following labels are",
				"not present in the model",
				"(use 'strict' = FALSE to ignore):",
				omxQuotes(missing))
			stop(msg)
		}
	}
	if (is.vector(lbound) && length(lbound) > 0 && all(sapply(lbound, is.na))) {
		lbound <- as.numeric(lbound)
	}
	if (is.vector(ubound) && length(ubound) > 0 && all(sapply(ubound, is.na))) {
		ubound <- as.numeric(ubound)
	}
	setParametersCheckVector(free, is.logical, 'free', 'logical')
	setParametersCheckVector(values, is.numeric, 'values', 'numeric')
	setParametersCheckVector(newlabels, is.character, 'newlabels', 'character')
	setParametersCheckVector(lbound, is.numeric, 'lbound', 'numeric')
	setParametersCheckVector(ubound, is.numeric, 'ubound', 'numeric')
	retval <- setParametersHelper(model, labels, free, values,
		newlabels, lbound, ubound, indep)
	if (!is.null(name)) {
		retval <- mxRename(retval, name)
	}
  retval@.modifiedSinceRun <- TRUE
	return(retval)
}

SBMatchHelper <- function(label, modelname) {
	if (!hasSquareBrackets(label)) {
		return(FALSE)
	}
	components <- splitSubstitution(label)
	fullname <- unlist(strsplit(components[[1]], imxSeparatorChar, fixed = TRUE))
	return(fullname[[1]] == modelname)
}

detectSBMatches <- function(model, labels) {
	modelname <- model@name
	targets <- which(sapply(labels, SBMatchHelper, modelname))
	return(targets)
}

setParametersHelper <- function(model, labels, free, values,
	newlabels, lbound, ubound, indep) {
	squarebrackets <- detectSBMatches(model, labels)
	model@matrices <- lapply(model@matrices, setParametersMatrix,
		labels, free, values, newlabels, lbound, ubound)
	if (length(squarebrackets) > 0) {
		model <- setSquareBracketsHelper(model, squarebrackets, labels, free, values, newlabels, lbound, ubound)
	}
	if(indep) {
		if (length(model@submodels) == 0) {
			return(model)
		}
		model@submodels <- lapply(model@submodels, setParametersHelper,
			labels, free, values, newlabels, lbound, ubound, indep)
	} else {
		select <- imxDependentModels(model)
		if (length(select) == 0) {
			return(model)
		}
		select <- lapply(select, setParametersHelper,
			labels, free, values, newlabels, lbound, ubound, indep)
		model@submodels <- c(select, imxIndependentModels(model))
	}
	return(model)
}

extractFirst <- function(x) {
	return(x[[1]])
}

extractSecond <- function(x) {
	return(x[[2]])
}

##' omxNameAnonymousParameters
##'
##' Assign new names to the unnamed parameters
##'
##' @param model the MxModel
##' @param indep whether models are independent
##' @return
##' a list with components for the new MxModel with named parameters, and the new names.
omxNameAnonymousParameters <- function(model, indep = FALSE) {
	rows <- lapply(model@matrices, getAnonymousRows)
    cols <- lapply(model@matrices, getAnonymousCols)
    newnames <- mapply(getAnonymousNames, rows)
	model@matrices <- mapply(assignAnonymousNames, model@matrices, rows, cols, newnames)
	newnames <- unlist(newnames)
	if(indep) {
		if (length(model@submodels) == 0) {
			return(list(model, newnames))
		}
		pairs <- lapply(model@submodels, omxNameAnonymousParameters, indep)
		submodels <- lapply(pairs, extractFirst)
		subnames <- unlist(lapply(pairs, extractSecond))
		names(submodels) <- names(model@submodels)
		model@submodels <- submodels
		newnames <- c(newnames, subnames)
	} else {
		select <- imxDependentModels(model)
		if (length(select) == 0) {
			return(list(model, newnames))
		}
		pairs <- lapply(select, omxNameAnonymousParameters, indep)
        submodels <- lapply(pairs, extractFirst)
        subnames <- unlist(lapply(pairs, extractSecond))
        names(submodels) <- names(select)
		model@submodels <- c(select, imxIndependentModels(model))
        newnames <- c(newnames, subnames)
	}
	return(list(model, newnames))
}

getAnonymousRows <- function(matrix) {
	select <- matrix@free & is.na(matrix@labels)
	if (imxSymmetricMatrix(matrix)) {
		select <- select & upper.tri(matrix@labels, diag=TRUE)
	}
	return(row(matrix@free)[select])
}

getAnonymousCols <- function(matrix) {
	select <- matrix@free & is.na(matrix@labels)
	if (imxSymmetricMatrix(matrix)) {
		select <- select & upper.tri(matrix@labels, diag=TRUE)
	}
	return(col(matrix@free)[select])
}

assignAnonymousNames <- function(matrix, rows, cols, newnames) {
	symmetry <- imxSymmetricMatrix(matrix)
	if (length(rows) > 0) {
		for(i in 1:length(rows)) {
			row <- rows[[i]]
			col <- cols[[i]]
			newname <- newnames[[i]]
			matrix@labels[row,col] <- newname
			if (symmetry) {
				matrix@labels[col,row] <- newname
			}
		}
	}
	return(matrix)
}

getAnonymousNames <- function(rows) {
	return(replicate(length(rows), imxUntitledName()))
}

omxAssignFirstParameters <- function(model, indep = FALSE) {
	params <- omxGetParameters(model, indep, fetch='all')
	if (!nrow(params)) return(model)
	model <- omxSetParameters(model, rownames(params),
		values = params$value, lbound = params$lbound, ubound = params$ubound, indep = indep)
  model
}

getParametersHelper <- function(matName, model, selection, fetch, labels) {
	amatrix <- model@matrices[[matName]]
	modelname <- model@name
	if (single.na(selection)) {
		select <- amatrix@free | !is.na(amatrix@labels)
	} else if (selection) {
		select <- amatrix@free
	} else {
		select <- !amatrix@free & !is.na(amatrix@labels)
	}
	if (length(labels)) for (lx in 1:length(labels)) {
		l1 <- labels[lx]
		if (hasSquareBrackets(l1)) {
			loc <- splitSubstitution(l1)
			fullname <- unlist(strsplit(loc[[1]], imxSeparatorChar, fixed = TRUE))
			if (fullname[[1]] == modelname && fullname[[2]] == matName) {
				select[as.integer(loc[2]), as.integer(loc[3])] <- TRUE
			}
		} else {
			select <- select | (!is.na(amatrix@labels) & amatrix@labels == l1)
		}
	}
	if (all(!select)) {
		return(numeric())
	}
	if (imxSymmetricMatrix(amatrix)) {
		triangle <- upper.tri(select, diag=TRUE)
		select <- select & triangle
	}
	theNames <- amatrix@labels[select]
	if (any(is.na(theNames))) {
		rows <- row(select)[select]
		cols <- col(select)[select]
		for(i in 1:length(theNames)) {
			if (is.na(theNames[[i]])) {
				theNames[[i]] <- paste(modelname, ".", amatrix@name,
				"[", rows[i], ",", cols[i], "]", sep ="")
			}
		}
	}
	if (fetch == "values") {
		theValues <- amatrix@values[select]
	} else if (fetch == "lbound") {
		theValues <- amatrix@lbound[select]
	} else if (fetch == "ubound") {
		theValues <- amatrix@ubound[select]
	} else if (fetch == "free") {
		theValues <- amatrix@free[select]
	}
	names(theValues) <- theNames
	return(theValues[!duplicated(theNames)])
}

#Need to use $ instead of @ accessor for all slots that can be condensed:
setParametersMatrix <- function(amatrix, names, free, values, newlabels, lbound, ubound) {
	labels <- amatrix@labels
	locations <- which(labels %in% names)
	indices <- match(labels[locations], names)
	if (!is.null(free)) {
		index2 <- ((indices - 1) %% length(free)) + 1
		amatrix$free[locations] <- as.logical(free[index2])
	}
	if (!is.null(values)) {
		index2 <- ((indices - 1) %% length(values)) + 1
		amatrix@values[locations] <- as.numeric(values[index2])
	}
	if (!is.null(newlabels)) {
		index2 <- ((indices - 1) %% length(newlabels)) + 1
		amatrix$labels[locations] <- as.character(newlabels[index2])
	}
	if (!is.null(lbound)) {
		index2 <- ((indices - 1) %% length(lbound)) + 1
		amatrix$lbound[locations] <- as.numeric(lbound[index2])
	}
	if (!is.null(ubound)) {
		index2 <- ((indices - 1) %% length(ubound)) + 1
		amatrix$ubound[locations] <- as.numeric(ubound[index2])
	}
	return(amatrix)
}

setSquareBracketsHelper <- function(model, squarebrackets, labels,
	free, values, newlabels, lbound, ubound) {
	for(i in 1:length(squarebrackets)) {
		nextbracket <- squarebrackets[[i]]
		nextlabel <- labels[[nextbracket]]
		components <- splitSubstitution(nextlabel)
		fullname <- unlist(strsplit(components[[1]], imxSeparatorChar, fixed = TRUE))
		matrixname <-fullname[[2]]
		row <- as.numeric(components[[2]])
		col <- as.numeric(components[[3]])
		amatrix <- model[[matrixname]]
		if (!is.null(amatrix) || !is(amatrix, "MxMatrix")) {
			isSymmetric <- imxSymmetricMatrix(amatrix)
			if (!is.null(free)) {
				index2 <- ((nextbracket - 1) %% length(free)) + 1
				amatrix$free[row,col] <- as.logical(free[index2])
				if (isSymmetric) {
					amatrix$free[col,row] <- as.logical(free[index2])
				}
			}
			if (!is.null(values)) {
				index2 <- ((nextbracket - 1) %% length(values)) + 1
				amatrix@values[row,col] <- as.numeric(values[index2])
				if (isSymmetric) {
					amatrix@values[col,row] <- as.numeric(values[index2])
				}
			}
			if (!is.null(newlabels)) {
				index2 <- ((nextbracket - 1) %% length(newlabels)) + 1
				amatrix$labels[row,col] <- as.character(newlabels[index2])
				if (isSymmetric) {
					amatrix$labels[col,row] <- as.character(newlabels[index2])
				}
			}
			if (!is.null(lbound)) {
				index2 <- ((nextbracket - 1) %% length(lbound)) + 1
				amatrix$lbound[row,col] <- as.numeric(lbound[index2])
				if (isSymmetric) {
					amatrix$lbound[col,row] <- as.numeric(lbound[index2])
				}
			}
			if (!is.null(ubound)) {
				index2 <- ((nextbracket - 1) %% length(ubound)) + 1
				amatrix$ubound[row,col] <- as.numeric(ubound[index2])
				if (isSymmetric) {
					amatrix$ubound[col,row] <- as.numeric(ubound[index2])
				}
			}
			model[[matrixname]] <- amatrix
		}
	}
	return(model)
}
