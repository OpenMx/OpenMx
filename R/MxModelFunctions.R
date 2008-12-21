omxGenerateMatrixList <- function(mxModel) {
	return(lapply(mxModel@matrices, generateMatrixListHelper))
}

omxGenerateSimpleMatrixList <- function(mxModel) {
	retval <- lapply(mxModel@matrices, generateMatrixListHelper)
	return(lapply(retval, as.matrix))
}

omxGenerateAlgebraList <- function(mxModel) {
    mList <- omxGenerateMatrixList(mxModel)
    retval <- lapply(mxModel@algebras, generateAlgebraHelper, names(mList))
    return(retval)
}

omxGenerateParameterList <- function(mxModel) {
	result <- list()
	if (length(mxModel@matrices) == 0) {
		return(result)
	}
	for(i in 1:length(mxModel@matrices)) {
		result <- generaterParameterListHelper(mxModel@matrices[[i]], result, i - 1)
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

omxUpdateModelValues <- function(mxModel, values) {
	pList <- omxGenerateParameterList(mxModel)
	if(length(pList) != length(values)) {
		stop(paste("This model has", length(pList), 
			"parameters, but you have given me", length(values),
			"values"))
	}
	if (length(pList) == 0) {
		return(mxModel)
	}
	for(i in 1:length(pList)) {
		mxModel <- updateModelValueHelper(pList[[i]], values[[i]], mxModel)
    }
	return(mxModel)
}

updateModelValueHelper <- function(triples, value, mxModel) {
	for(i in 1:length(triples)) {
		triple <- triples[[i]]
		mat <- triple[1] + 1
		row <- triple[2]
		col <- triple[3]
		mxModel@matrices[[mat]]@values[row,col] <- value			
	}
	return(mxModel)
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
