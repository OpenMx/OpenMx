omxGenerateMatrixList <- function(mxModel) {
	return(lapply(mxModel@matrices, generateMatrixListHelper))
}

omxGenerateSimpleMatrixList <- function(mxModel) {
	retval <- lapply(mxModel@matrices, generateMatrixListHelper)
	return(lapply(retval, as.matrix))
}

omxGenerateAlgebraList <- function(mxModel) {
    retval <- lapply(mxModel@algebras, generateAlgebraHelper, 
    	names(mxModel@matrices), names(mxModel@algebras))
    return(retval)
}

omxGenerateParameterList <- function(mxModel) {
	result <- list()
	if (length(mxModel@matrices) == 0) {
		return(result)
	}
	for(i in 1:length(mxModel@matrices)) {
		result <- generaterParameterListHelper(
			mxModel@matrices[[i]], 
			result, i - 1)
	}	
	return(result)
}

omxGenerateValueList <- function(mxModel) {
	mList <- omxGenerateMatrixList(mxModel)
	pList <- omxGenerateParameterList(mxModel)
	retval <- vector()
	if (length(pList) == 0) {
		return(retval)
	}
	for(i in 1:length(pList)) {
		parameter <- pList[[i]]
		if (length(parameter) > 1) {
			values <- sapply(parameter, generateValueHelper, mList)
			if (!all(values == values[[1]])) {
				warning(paste('Parameter',i,'has multiple start values.',
					'Selecting', values[[1]]))
			}
			retval[i] <- values[[1]]
		} else {
			retval[i] <- generateValueHelper(parameter[[1]], mList)
		}
    }
	return(retval)	
}

generateValueHelper <- function(triple, mList) {
	mat <- triple[1] + 1
	row <- triple[2]
	col <- triple[3]
	return(mList[[mat]][row,col])
}

omxUpdateModelValues <- function(treeModel, flatModel, values) {
	pList <- omxGenerateParameterList(flatModel)
	if(length(pList) != length(values)) {
		stop(paste("This model has", length(pList), 
			"parameters, but you have given me", length(values),
			"values"))
	}
	if (length(pList) == 0) {
		return(treeModel)
	}
	for(i in 1:length(pList)) {
		treeModel <- updateModelValueHelper(
			pList[[i]], values[[i]], treeModel, flatModel)
    }
	return(treeModel)
}

updateModelValueHelper <- function(triples, value, treeModel, flatModel) {
	for(i in 1:length(triples)) {
		triple <- triples[[i]]
		mat <- triple[1] + 1
		row <- triple[2]
		col <- triple[3]
		name <- flatModel@matrices[[mat]]@name
		if(!is.null(treeModel[[name]])) {
			treeModel[[name]]@values[row,col] <- value
		}
	}
	treeModel@submodels <- lapply(treeModel@submodels, 
		function(x) { updateModelValueHelper(
			triples, value, x, flatModel) })
	return(treeModel)
}

omxUpdateModelAlgebras <- function(treeModel, flatModel, values) {
	aList <- omxGenerateAlgebraList(flatModel)
	if(length(aList) != length(values)) {
		stop(paste("This model has", length(aList), 
			"algebras, but you have given me", length(values),
			"values"))
	}
	if (length(aList) == 0) {
		return(treeModel)
	}	
	treeModel <- updateModelAlgebraHelper(aList, values, treeModel)
	return(treeModel)
}

updateModelAlgebraHelper <- function(aList, values, model) {
	aNames <- names(aList)
	for(i in 1:length(aList)) {
		name <- aNames[[i]]
		candidate <- model[[name]]
		if (!is.null(candidate) && is(candidate,"MxAlgebra")) {
			model[[name]]@result <- Matrix(values[[i]])
		}
	}
	model@submodels <- lapply(model@submodels, function(x) {updateModelAlgebraHelper(aList, values, x)})
	return(model)
}


omxUpdateModelObjective <- function(model, result) {
	if(!is.null(model@objective)) {
		model@objective@result <- Matrix(result)
	}
	return(model)
}

omxLocateIndex <- function(model, name) {
	matrixNumber <- match(name, names(model@matrices))
	algebraNumber <- match(name, names(model@algebras))
	if (is.na(matrixNumber) && is.na(algebraNumber)) {
		return(NA)
	} else if (is.na(algebraNumber)) {
		return(- matrixNumber)
	} else {
		return(algebraNumber - 1)
	}
}
