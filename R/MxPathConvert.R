convertModelA <- function(model) {
	variables <- c(model@manifestVars, model@latentVars)
	len <- length(variables)
	names <- list(variables, variables)
	values <- Matrix(0, nrow = len, ncol = len, dimnames = names)
	specification <- new("MxSparseMatrix", 
		nrow = len, ncol = len, dimnames = names)
	npaths <- dim(model@paths)[[1]]
	for(i in 1:npaths) {
		apath <- model@paths[i,]
		values[apath['to'][[1]], apath['from'][[1]]] <- getValuesA(apath)
		specification[apath['to'][[1]], apath['from'][[1]]] <- getSpecificationA(apath)
	}
	retval <- mxMatrix("Full", values, specification, nrow = len, ncol = len)
	return(retval)
}

getValuesA <- function(apath) {
	if (apath[['arrows']] == 1) {
 		if (is.null(apath[['startVal']])) {
			startVal <- 1.00
		} else {
			startVal <- apath[['startVal']]
		}
		return(startVal)
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
	specification <- new("MxSymmetricSparse", nrow = len, ncol = len, dimnames = names)
	npaths <- dim(model@paths)[[1]]
	for(i in 1:npaths) {
		apath <- model@paths[i,]
		values[apath['to'][[1]], apath['from'][[1]]] <- getValuesS(apath)
		specification[apath['to'][[1]], apath['from'][[1]]] <- getSpecificationS(apath)
	}
	retval <- mxMatrix("Symm", values, specification, nrow = len, ncol = len)
	return(retval)
}

getValuesS <- function(apath) {
	if (apath[['arrows']] == 2) {
 		if (is.null(apath[['startVal']])) {
			startVal <- 1.00
		} else {
			startVal <- apath[['startVal']]
		}
		return(startVal)
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