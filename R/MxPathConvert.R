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
	values <- Matrix(0, nrow = len, ncol = len, dimnames = names)
	spec <- new("MxSparseMatrix", nrow = len, ncol = len, 
		dimnames = names)
	npaths <- dim(model@paths)[[1]]
	for(i in 1:npaths) {
		apath <- model@paths[i,]
		values[apath['to'][[1]], apath['from'][[1]]] <- getValuesA(apath)
		spec[apath['to'][[1]], apath['from'][[1]]] <- getSpecificationA(apath)
	}
	retval <- mxMatrix("Full", values, spec)
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

getSpecificationA <- function(apath) {
	if ((apath[['arrows']] == 1) && (apath[['free']] == TRUE)) {
		return(NA)
	} else {
		return(0)
	}	
}


convertModelS <- function(model) {
	variables <- c(model@manifestVars, model@latentVars)
	len <- length(variables)
	names <- list(variables, variables)
	values <- new("MxSymmetricSparse", nrow = len, ncol = len, dimnames = names)
	spec <- new("MxSymmetricSparse", nrow = len, ncol = len, dimnames = names)
	npaths <- dim(model@paths)[[1]]
	for(i in 1:npaths) {
		apath <- model@paths[i,]
		values[apath['to'][[1]], apath['from'][[1]]] <- getValuesS(apath)
		spec[apath['to'][[1]], apath['from'][[1]]] <- getSpecificationS(apath)
	}
	retval <- mxMatrix("Symm", values, spec, nrow = len, ncol = len)
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

getSpecificationS <- function(apath) {
	if ((apath[['arrows']] == 2) && (apath[['free']] == TRUE)) {
		return(NA)
	} else {
		return(0)
	}	
}

convertModelF <- function(model) {
	variables <- c(model@manifestVars, model@latentVars)
	len <- length(variables)
	names <- list(model@manifestVars, variables)
	matValues <- diag(nrow = length(model@manifestVars), ncol = len)
	values <- Matrix(matValues, dimnames = names)
	retval <- mxMatrix("Full", values)
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
