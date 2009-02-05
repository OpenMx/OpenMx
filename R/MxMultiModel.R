omxCheckNamespace <- function(model) {
	omxCheckNamespaceHelper(model, c(model@name))
}

omxNamespaceErrorMessage <- function(rlist) {
	if (length(rlist) == 1) {
		return(paste("The name", omxQuotes(rlist), 
		"appears more than once in the model\n"))
	} else {
		return(paste("The names", omxQuotes(rlist), 
		"appear more than once in the model\n"))
	}
}

omxGetNames <- function(lst) {
	lapply(lst, function(x) { x@name } )	
}

omxNameAlignment <- function(lst1, lst2) {
	return(paste("The names", omxQuotes(lst1[lst1 != lst2]),
		"do not match their designations"))
}

omxCheckNamespaceHelper <- function(model, nlist) {
	omxNameAlignment(names(model@matrices), omxGetNames(model@matrices))
	matrixIntersect <- intersect(names(model@matrices), nlist)
	if (length(matrixIntersect) > 0) {
		stop(omxNamespaceErrorMessage(matrixIntersect), call.=FALSE)
	}
	nlist <- append(nlist, names(model@matrices))
	omxNameAlignment(names(model@algebras), omxGetNames(model@algebras))
	algebraIntersect <- intersect(names(model@algebras), nlist)
	if (length(algebraIntersect) > 0) {
		stop(omxNamespaceErrorMessage(algebraIntersect), call.=FALSE)
	}
	nlist <- append(nlist, names(model@algebras))
	omxNameAlignment(names(model@submodels), omxGetNames(model@submodels))
	modelIntersect <- intersect(names(model@submodels), nlist)
	if (length(modelIntersect) > 0) {
		stop(omxNamespaceErrorMessage(modelIntersect), call.=FALSE)
	}
	nlist <- append(nlist, names(model@submodels))
	if (!is.null(model@objective) && (model@objective@name %in% nlist)) {
		stop(omxNamespaceErrorMessage(model@objective@name), call.=FALSE)
	} else if(!is.null(model@objective)) {
		nlist <- append(nlist, model@objective@name)
	}
	if (length(model@submodels) > 0) {
		for(i in 1:length(model@submodels)) {
			nlist <- omxCheckNamespaceHelper(model@submodels[[i]], nlist)
		}
	}
	return(nlist)
}

omxShareData <- function(model) {
	return(omxShareDataHelper(model, model@data))
}

omxShareDataHelper <- function(model, current) {
	if((length(model@data) == 0) && (model@independent == TRUE)) {
		model@data <- current
	} else {
		current <- model@data
	}
	submodels <- lapply(model@submodels, function(x) { omxShareDataHelper(x, current) })
	model@submodels <- submodels
	return(model)
}

omxGetIndependents <- function(model) {
	return(omxGetIndependentsHelper(model, list()))
}

omxGetIndependentsHelper <- function(model, lst) {
	indep <- lapply(model@submodels, function(x) {
		if(x@independent == TRUE) {x}
		else {NULL}
	})
	res <- append(lst, indep[!sapply(indep, is.null)])
	if (length(model@submodels) > 0) {
		for(i in 1:length(model@submodels)) {
			res <- omxGetIndependentsHelper(model@submodels[[i]], res)
		}
	}
	return(res)
}

omxFreezeMatrix <- function(mxMatrix) {
	rows <- nrow(mxMatrix)
	cols <- ncol(mxMatrix)
	if (is(mxMatrix@specification, "MxSymmetricSparse")) {
		mxMatrix@specification <- new("MxSymmetricSparse", 
			nrow = rows, ncol = cols)
	} else {
		mxMatrix@specification <- new("MxSparseMatrix", 
			nrow = rows, ncol = cols)
	}
	return(mxMatrix)
}

omxFreezeAlgebra <- function(mxAlgebra) {
	if(is.null(mxAlgebra@result)) return(NULL)
	res <- mxMatrix(values = mxAlgebra@result, 
		name = mxAlgebra@name)
	return(res)
}

omxFreezeObjective <- function(model) {
	objective <- model@objective
	if (!is.null(objective)) {
		model[[objective@name]] <- NULL
		if (!is.null(objective@result)) {
			model[[objective@name]] <- mxMatrix(
				values = objective@result,
				name = objective@name)
		}
	}
	return(model)
}

omxFreezeModel <- function(model) {
	model <- omxFreezeObjective(model)
	model@matrices <- lapply(model@matrices, omxFreezeMatrix)
	algebras <- lapply(model@algebras, omxFreezeAlgebra)
	algebras <- algebras[!sapply(algebras, is.null)]
	model@matrices <- append(model@matrices, algebras)
	model@algebras <- list()
	model@submodels <- lapply(model@submodels, omxFreezeModel)
	return(model)
}