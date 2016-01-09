#
#   Copyright 2007-2016 The OpenMx Project
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
	# add '+' on both sides to ensure that we always get 2 pieces
	split2 <- unlist(strsplit(paste("+",split1[[2]],"+",sep=""), ",", fixed = TRUE))
	row <- substr(split2[1], 2, nchar(split2[1]))
	col <- substr(split2[2], 1, nchar(split2[2])-1)
	return(c(identifier, row, col))
}

findSquareBrackets <- function(aMatrix) {
	aMatrix@.squareBrackets <- apply(aMatrix@labels, c(1,2), hasSquareBrackets)
	return(aMatrix)
}

hasSquareBrackets <- function(input) {
	if (is.na(input)) {
		return(FALSE)
	}
    match1 <- grep("[", input, fixed=TRUE)
	if (length(match1) == 0) {
		return(FALSE)
	}
    match2 <- grep("]", input, fixed=TRUE)
    return(length(match2) > 0)
}

verifySquareBracket <- function(squareBracket, matrixName) {
	components <- splitSubstitution(squareBracket)
	row <- components[[2]]
	col <- components[[3]]
	rowmatch <- grep("^[0-9]+$", row)
	colmatch <- grep("^[0-9]+$", col)
	if (length(rowmatch) == 0 || length(colmatch) == 0) {
		msg <- paste("Illegal label",
			omxQuotes(squareBracket),
			"detected in matrix", paste(omxQuotes(matrixName), '.', sep=''),
			"Square brackets must contain numeric literals",
			"when used inside of labels.")
		stop(msg, call. = FALSE)
	}
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
		select <- matrix@.squareBrackets
		rows <- row(labels)[select]
		cols <- col(labels)[select]
		subs <- labels[select]
		if (length(subs) > 0) {
			for (j in 1:length(subs)) {
				components <- splitSubstitution(subs[[j]])
				identifier <- components[[1]]
				fromrow <- as.integer(components[[2]]) - 1L
				fromcol <- as.integer(components[[3]]) - 1L
				torow <- as.integer(rows[j]) - 1L
				tocol <- as.integer(cols[j]) - 1L
				index <- imxLocateIndex(model, identifier, name)
				len <- length(retval[[name]])
				retval[[name]][[len + 1]] <- c(index, fromrow, fromcol, torow, tocol)
			}
		}
	}
	return(retval)
}
