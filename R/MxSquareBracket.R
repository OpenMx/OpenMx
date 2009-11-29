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

splitSubstitution <- function(input) {
	split1 <- unlist(strsplit(input, "[\\[\\]]", perl = TRUE))
	identifier <- split1[[1]]
	split2 <- unlist(strsplit(split1[[2]], ",", fixed = TRUE))
	row <- split2[[1]]
	col <- split2[[2]]
	return(c(identifier, row, col))
}

isSubstitution <- function(input) {
    match <- grep("^([^\\[\\]])+\\[[0-9]+,[0-9]+\\]$", input, perl = TRUE, value = TRUE)
    return(length(match) > 0)
}

hasSquareBrackets <- function(input) {
    match <- grep("[\\[\\]]", input, perl = TRUE, value = TRUE)
    return(length(match) > 0)
}

checkSquareBracketEvaluation <- function(model, flatModel, labelsData) {
	if(length(flatModel@matrices) == 0) { return() }
	for(i in 1:length(flatModel@matrices)) {
		checkSquareBracketMatrix(flatModel@matrices[[i]], model, flatModel, labelsData)
	}
}

checkSquareBracketMatrix <- function(matrix, model, flatModel, labelsData) {
	labels <- matrix@labels
	labels <- labels[!is.na(labels)]
	if (length(labels) == 0) { return() }
	subs <- sapply(labels, isSubstitution)
	labels <- labels[subs]
	if (length(labels) == 0) { return() }
	for(i in 1:length(labels)) {
		substitute <- labels[[i]]
		pieces <- splitSubstitution(substitute)
		identifier <- pieces[[1]]
		idenrow <- as.numeric(pieces[[2]])
		idencol <- as.numeric(pieces[[3]])
		entity <- flatModel[[identifier]]
		if (is.null(entity)) {
			stop(paste("Unknown reference", 
 				omxQuotes(simplifyName(identifier, model@name)),
				"detected in the matrix", omxQuotes(simplifyName(matrix@name, model@name)),
				"in model", omxQuotes(model@name)), call. = FALSE)
		}
		if (is(entity, "MxMatrix")) {
			subrow <- nrow(entity)
			subcol <- ncol(entity)
		} else if (is(entity, "MxAlgebra")) {
			value <- as.matrix(eval(computeSymbol(as.symbol(identifier), flatModel, labelsData)))
			subrow <- nrow(value)
			subcol <- ncol(value)
		} else if (is(entity, "MxObjectiveFunction")) {
			subrow <- 1
			subcol <- 1
		} else {
			stop(paste("Cannot apply the substitution using", 
 				omxQuotes(simplifyName(identifier, model@name)),
				"detected in the matrix", omxQuotes(simplifyName(matrix@name, model@name)),
				"in model", omxQuotes(model@name)), call. = FALSE)			
		}
		if (idenrow < 0 || idencol < 0 || idenrow > subrow || idencol > subcol) {
			identifier <- simplifyName(identifier, model@name)
			substitute <- paste(identifier, '[', idenrow, ',', idencol, ']', sep = '')
			stop(paste("The substitution", 
 				omxQuotes(substitute),
				"detected in the matrix", omxQuotes(simplifyName(matrix@name, model@name)),
				"in model", omxQuotes(model@name),
				"has invalid (row,col) values"), call. = FALSE)
		}
	}
}

omxComputeSubstitution <- function(matrix, model) {
	labels <- matrix@labels
	select <- !apply(labels, c(1,2), is.na) & apply(labels, c(1,2), isSubstitution)
	rows <- row(labels)[select]
	cols <- col(labels)[select]
	subs <- labels[select]
	if (length(subs) == 0) { return(matrix@values) }
	value <- matrix@values
	for(i in 1:length(subs)) {
		pieces <- splitSubstitution(subs[[i]])
		identifier <- pieces[[1]]
		sourcerow <- as.numeric(pieces[[2]])
		sourcecol <- as.numeric(pieces[[3]])
		other <- eval(substitute(mxEval(x, model, compute=TRUE), list(x = as.symbol(identifier)))) 
		value[rows[[i]],cols[[i]]] <- other[sourcerow,sourcecol]
	}
	return(value)
}

generateMatrixReferences <- function(model) {
	matnames <- names(model@matrices)
	retval <- replicate(length(matnames), list())
	names(retval) <- matnames
	if (length(model@matrices) == 0) {
		return(retval)
	}
	for (i in 1:length(model@matrices)) {
		matrix <- model@matrices[[i]]
		name <- matrix@name
		labels <- matrix@labels
		select <- !apply(labels, c(1,2), is.na) & apply(labels, c(1,2), isSubstitution)
		rows <- row(labels)[select]
		cols <- col(labels)[select]
		subs <- labels[select]
		if (length(subs) > 0) {
			for (j in 1:length(subs)) {
				pieces <- splitSubstitution(subs[[j]])
				identifier <- pieces[[1]]
				fromrow <- as.integer(pieces[[2]]) - 1L
				fromcol <- as.integer(pieces[[3]]) - 1L
				torow <- as.integer(rows[j] - 1)
				tocol <- as.integer(cols[j] - 1)
				index <- omxLocateIndex(model, identifier, name)
				len <- length(retval[[name]])
				retval[[name]][[len + 1]] <- c(index, fromrow, fromcol, torow, tocol)
			}
		}
	}
	return(retval)
}
