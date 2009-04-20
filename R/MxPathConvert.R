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


convertModelA <- function(model) {
	variables <- c(model@manifestVars, model@latentVars)
	len <- length(variables)
	names <- list(variables, variables)
	values <- matrix(0, len, len, dimnames = names)
	free <- matrix(FALSE, len, len, dimnames = names)
	labels <- matrix(as.character(NA), len, len, dimnames = names)
	npaths <- dim(model@paths)[[1]]
	for(i in 1:npaths) {
		apath <- model@paths[i,]
		values[apath['to'][[1]], apath['from'][[1]]] <- getValuesA(apath)
		free[apath['to'][[1]], apath['from'][[1]]] <- getFreeA(apath)
		labels[apath['to'][[1]], apath['from'][[1]]] <- getLabelsA(apath)
	}
	retval <- mxMatrix("Full", values, free, labels, name ="A")
	return(retval)
}

getValuesA <- function(apath) {
	if (apath[['arrows']] == 1) {
 		if (is.null(apath[['start']])) {
			start <- 1.00
		} else {
			start <- apath[['start']]
		}
		return(start)
	} else {
		return(0.0)
	}
}

getFreeA <- function(apath) {
	if ((apath[['arrows']] == 1) && (apath[['free']] == TRUE)) {
		return(TRUE)
	} else {
		return(FALSE)
	}	
}

getLabelsA <- function(apath) {
	if ((apath[['arrows']] == 1) && (apath[['free']] == TRUE)) {
		if (is.null(apath[['name']])) {
			return(as.character(NA))
		} else {
			return(apath[['name']])
		}
	} else {
		return(as.character(NA))
	}
}


convertModelS <- function(model) {
	variables <- c(model@manifestVars, model@latentVars)
	len <- length(variables)
	names <- list(variables, variables)
	values <- matrix(0, len, len, dimnames = names)	
	free <- matrix(FALSE, len, len, dimnames = names)
	labels <- matrix(as.character(NA), len, len, dimnames = names)
	npaths <- dim(model@paths)[[1]]
	for(i in 1:npaths) {
		apath <- model@paths[i,]
		pathValue <- getValuesS(apath)
		pathFree <- getFreeS(apath)
		pathLabel <- getLabelsS(apath)
		values[apath['to'][[1]], apath['from'][[1]]] <- pathValue
		values[apath['from'][[1]], apath['to'][[1]]] <- pathValue
		free[apath['to'][[1]], apath['from'][[1]]] <- pathFree
		free[apath['from'][[1]], apath['to'][[1]]] <- pathFree
		labels[apath['to'][[1]], apath['from'][[1]]] <- pathLabel
		labels[apath['from'][[1]], apath['to'][[1]]] <- pathLabel
	}
	retval <- mxMatrix("Symm", values, free, labels, name = "S")
	return(retval)
}

getValuesS <- function(apath) {
	if (apath[['arrows']] == 2) {
 		if (is.null(apath[['start']])) {
			start <- 1.00
		} else {
			start <- apath[['start']]
		}
		return(start)
	} else {
		return(0.0)
	}
}

getFreeS <- function(apath) {
	if ((apath[['arrows']] == 2) && (apath[['free']] == TRUE)) {
		return(TRUE)
	} else {
		return(FALSE)
	}	
}

getLabelsS <- function(apath) {
	if ((apath[['arrows']] == 2) && (apath[['free']] == TRUE)) {
		if (is.null(apath[['name']])) {
			return(as.character(NA))
		} else {
			return(apath[['name']])
		}
	} else {
		return(as.character(NA))
	}	
}


convertModelF <- function(model) {
	variables <- c(model@manifestVars, model@latentVars)
	len <- length(variables)
	names <- list(model@manifestVars, variables)
	values <- diag(nrow = length(model@manifestVars), ncol = len)
	dimnames(values) <- names
	retval <- mxMatrix("Full", values, name = "F")
	return(retval)
}

omxConvertPathModel <- function(model) {
	if(!inherits(try(omxLocateIndex(model, "A", model@name), 
		silent = TRUE), "try-error")) {
		warning("Overwriting existing 'A' matrix/algebra while converting RAM-style model")
	}
	if(!inherits(try(omxLocateIndex(model, "S", model@name),
		silent = TRUE), "try-error")) {
		warning("Overwriting existing 'S' matrix/algebra while converting RAM-style model")
	}
	if(!inherits(try(omxLocateIndex(model, "F", model@name),
		silent = TRUE), "try-error")) {
		warning("Overwriting existing 'F' matrix/algebra while converting RAM-style model")
	}
	model[['A']] <- convertModelA(model)
	model[['S']] <- convertModelS(model)
	model[['F']] <- convertModelF(model)
	return(model)
}
