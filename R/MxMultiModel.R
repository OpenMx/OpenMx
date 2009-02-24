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

omxCheckNamedEntity <- function(model, slotname, nlist) {
	entity <- slot(model, slotname)
	omxNameAlignment(names(entity), omxGetNames(entity))
	entityIntersect <- intersect(names(entity), nlist)
	if (length(entityIntersect) > 0) {
		stop(omxNamespaceErrorMessage(entityIntersect), call.=FALSE)
	}
	nlist <- append(nlist, names(entity))
	return(nlist)
}

omxCheckNamespaceHelper <- function(model, nlist) {
	nlist <- omxCheckNamedEntity(model, "matrices", nlist)
	nlist <- omxCheckNamedEntity(model, "algebras", nlist)
	nlist <- omxCheckNamedEntity(model, "submodels", nlist)
	nlist <- omxCheckNamedEntity(model, "constraints", nlist)
	nlist <- omxCheckNamedEntity(model, "bounds", nlist)
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
	if(is.null(model@data) && (model@independent == TRUE)) {
		model@data <- current
	} else if (!is.null(model@data)) {
		current <- model@data
	}
	submodels <- lapply(model@submodels, function(x)
		{ omxShareDataHelper(x, current) })
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

omxFlattenModel <- function(model) {
	res <- new("MxFlatModel", model, list(), list())
	if (!is.null(model@objective)) {
		defaultData <- model@objective@name
	} else {
		defaultData <- model@name
		res@datasets[[defaultData]] <- model@data
	}
	res <- omxFlattenModelHelper(model, res, defaultData)
	res@submodels <- list()
	return(res)
}

omxFlattenModelHelper <- function(model, dest, defaultData) {
	if (!is.null(model@objective)) {
		name <- model@objective@name
		dest@objectives[[name]] <- model@objective
		if(is.null(model@data)) {
			dest@datasets[[name]] <- defaultData
		} else {
			dest@datasets[[name]] <- model@data
		}
	}
	if (is.null(dest@datasets[[defaultData]]) && !is.null(model@data)) {
		if (!is.null(model@objective)) {
			defaultData <- model@objective@name
		} else {
			defaultData <- model@name
			dest@datasets[[defaultData]] <- model@data
		}
	}
	if (length(model@submodels) > 0) {
		for(i in 1:length(model@submodels)) {
			submodel <- model@submodels[[i]]
			dest@matrices    <- append(dest@matrices, submodel@matrices)
			dest@algebras    <- append(dest@algebras, submodel@algebras)
			dest@constraints <- append(dest@constraints, 
				submodel@constraints) 
			dest@bounds      <- append(dest@bounds, submodel@bounds)
			dest <- omxFlattenModelHelper(submodel, dest, defaultData)
		}
	}
	return(dest)
}

omxReplaceModels <- function(model, replacements) {
	mnames <- names(model@submodels)
	rnames <- names(replacements)
	inames <- intersect(mnames, rnames)
	if (length(inames) > 0) {
		for(i in 1:length(inames)) {
			name <- inames[[i]]
			model[[name]] <- replacements[[name]]
		}
	}
	if (length(model@submodels) > 0) {
		for(i in 1:length(model@submodels)) {
			model@submodels[[i]] <- 
				omxReplaceModels(model@submodels[[i]], replacements)
		}
	}
	return(model)
}