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

generateMatrixReferences <- function(model) {
	matnames  <- names(model@matrices)
	retval <- replicate(length(matnames), list())
	names(retval) <- matnames
	if (length(model@matrices) == 0) {
		return(retval)
	}
	for (i in 1:length(model@matrices)) {
		matrix <- model@matrices[[i]]
		name <- matrix@name
		labels <- matrix@labels
		select <- !is.na(labels)
		if (any(select == TRUE)) {
			parameterNames <- labels[select]
			rows <- row(labels)[select]
			cols <- col(labels)[select]			
			for (j in 1:length(parameterNames)) {
				parameterName <- parameterNames[j]
				row <- rows[j] - 1
				col <- cols[j] - 1
				if (parameterName %in% matnames) {
					len <- length(retval[[parameterName]])
					retval[[parameterName]][[len + 1]] <- c(i - 1, row, col)
				}
			}
		}
	}
	return(retval)
}

generateAlgebraReferences <- function(model, algebras, objectives) {
	values <- c(algebras, objectives)
	references <- generateAlgebraReferencesHelper(model)
	retval <- mapply(function(x,y) { c(list(x), y) }, 
		values, references, SIMPLIFY = FALSE)
	return(retval)
}

generateAlgebraReferencesHelper <- function(model) {
	matnames  <- names(model@matrices)
	algnames  <- c(names(model@algebras), names(model@objectives)) 
	retval <- replicate(length(algnames), list())
	names(retval) <- algnames
	if (length(model@matrices) == 0) {
		return(retval)
	}
	for (i in 1:length(model@matrices)) {
		matrix <- model@matrices[[i]]
		name <- matrix@name
		labels <- matrix@labels
		select <- !is.na(labels)
		if (any(select == TRUE)) {
			parameterNames <- labels[select]
			rows <- row(labels)[select]
			cols <- col(labels)[select]			
			for (j in 1:length(parameterNames)) {
				parameterName <- parameterNames[j]
				row <- rows[j] - 1
				col <- cols[j] - 1
				if (parameterName %in% algnames) {
					len <- length(retval[[parameterName]])
					retval[[parameterName]][[len + 1]] <- c(i - 1, row, col)
				}
			}
		}
	}
	return(retval)
}
