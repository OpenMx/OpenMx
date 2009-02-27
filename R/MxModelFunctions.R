omxGenerateMatrixList <- function(mxModel) {
	return(lapply(mxModel@matrices, generateMatrixListHelper))
}

omxGenerateSimpleMatrixList <- function(mxModel) {
	retval <- lapply(mxModel@matrices, generateMatrixListHelper)
	return(lapply(retval, as.matrix))
}


omxGenerateAlgebraList <- function(mxModel) {
	mNames <- names(mxModel@matrices)
	aNames <- names(mxModel@algebras)
	oNames <- names(mxModel@objectives)
    retval <- lapply(mxModel@algebras, generateAlgebraHelper, 
    	mNames, append(aNames, oNames))
    return(retval)
}

omxGenerateParameterList <- function(mxModel) {
	result <- list()
	if (length(mxModel@matrices) == 0) {
		return(result)
	}
	for(i in 1:length(mxModel@matrices)) {
		result <- omxGenerateParameterListHelper(
			mxModel@matrices[[i]], 
			mxModel@bounds, result, i - 1)
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
		parameter <- parameter[3:length(parameter)] # Remove (min, max) bounds
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

omxConvertObjectives <- function(flatModel) {
	retval <- lapply(flatModel@objectives, function(x) {
		omxObjFunConvert(x, flatModel)
	})
	return(retval)
}

omxObjectiveIndex <- function(flatModel) {
	objective <- flatModel@objective
	if(is.null(objective)) {
		return(NULL)
	} else {
		return(omxLocateIndex(flatModel, objective@name, flatModel@name))
	}
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
		parameters <- pList[[i]]
		parameters <- parameters[3:length(parameters)] # Remove (min, max) bounds
		treeModel <- updateModelValueHelper(
			parameters, values[[i]], treeModel, flatModel)
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
	aNames <- names(flatModel@algebras)
	oNames <- names(flatModel@objectives)
	aList <- append(aNames, oNames)
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
	print(length(aList))
	print(values)
	for(i in 1:length(aList)) {
		name <- aList[[i]]
		candidate <- model[[name]]
		if (!is.null(candidate) && !is.nan(values[[i]]) && 
			(is(candidate,"MxAlgebra") || (is(candidate,"MxObjective")))) {
			model[[name]]@result <- Matrix(values[[i]])
		}
	}
	model@submodels <- lapply(model@submodels, function(x) {
		updateModelAlgebraHelper(aList, values, x)})
	return(model)
}

omxLocateIndex <- function(model, name, referant) {
	mNames <- names(model@matrices)
	aNames <- names(model@algebras)
	oNames <- names(model@objectives)		
	matrixNumber <- match(name, mNames)
	algebraNumber <- match(name, append(aNames, oNames))
	if (is.na(matrixNumber) && is.na(algebraNumber)) {
		msg <- paste("The reference", omxQuotes(name),
			"does not exist.  It is used by the named entity",
			omxQuotes(referant),".")
		stop(msg, call.=FALSE)
	} else if (is.na(algebraNumber)) {
		return(- matrixNumber)
	} else {
		return(algebraNumber - 1)
	}
}
