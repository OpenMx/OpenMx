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

hasSquareBrackets <- function(input) {
    match <- grep("[\\[\\]]", input, perl = TRUE, value = TRUE)
    return(length(match) > 0)
}

convertSquareBracketLabels <- function(model) {
	if(length(model@matrices) > 0) {
		for(i in 1:length(model@matrices)) {
			model <- convertSquareBracketHelper(model, i)
		}
	}
	if(length(model@submodels) > 0) {
		for(i in 1:length(model@submodels)) {
			if (model@submodels[[i]]@independent == FALSE) {
				model@submodels[[i]] <- convertSquareBracketLabels(model@submodels[[i]])
			}
		}
	}
	return(model)
}

undoSquareBracketLabels <- function(model) {
	if(length(model@matrices) > 0) {
		for(i in 1:length(model@matrices)) {
			model <- undoSquareBracketHelper(model, i)
		}
	}
	if(length(model@submodels) > 0) {
		for(i in 1:length(model@submodels)) {
			if (model@submodels[[i]]@independent == FALSE) {
				model@submodels[[i]] <- undoSquareBracketLabels(model@submodels[[i]])
			}
		}
	}
	return(model)
}

convertSquareBracketHelper <- function(model, index) {
	target <- model@matrices[[index]]
	labels <- target@labels
	select <- !apply(labels, c(1,2), is.na) & apply(labels, c(1,2), hasSquareBrackets)
	rows <- row(labels)[select]
	cols <- col(labels)[select]
	subs <- labels[select]
	if (length(subs) > 0) {
		for (i in 1:length(subs)) {
			row <- rows[[i]]
			col <- cols[[i]]
			formula <- parse(text = subs[[i]])
			name <- paste(target@name, '[', row, ',', col, ']', sep='')
			algebra <- eval(substitute(mxAlgebra(x), list(x = formula[[1]])))
			model[[name]] <- algebra
			model@matrices[[index]]@labels[row,col] <- name
		}
	}
	return(model)
}

undoSquareBracketHelper <- function(model, index) {
	target <- model@matrices[[index]]
	labels <- target@labels
	select <- !apply(labels, c(1,2), is.na) & apply(labels, c(1,2), hasSquareBrackets)
	rows <- row(labels)[select]
	cols <- col(labels)[select]
	subs <- labels[select]
	if (length(subs) > 0) {
		for (i in 1:length(subs)) {
			row <- rows[[i]]
			col <- cols[[i]]
			name <- paste(target@name, '[', row, ',', col, ']', sep='')
			algebra <- model[[name]]
			newlabel <- deparse(algebra@formula, width.cutoff = 500)
			model@matrices[[index]]@labels[row,col] <- newlabel
			model[[name]] <- NULL
		}
	}
	return(model)
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
		select <- !apply(labels, c(1,2), is.na) & apply(labels, c(1,2), hasSquareBrackets)
		rows <- row(labels)[select]
		cols <- col(labels)[select]
		subs <- labels[select]
		if (length(subs) > 0) {
			for (j in 1:length(subs)) {
				identifier <- subs[[j]]
				fromrow <- 0L
				fromcol <- 0L
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
