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


generateMatrixList <- function(mxModel) {
	matvalues <- lapply(mxModel@matrices, generateMatrixValuesHelper)
	matnames  <- names(mxModel@matrices)
	names(matvalues) <- matnames
	references <- generateMatrixReferences(mxModel)
	retval <- mapply(function(x,y) { c(list(x), y) }, 
			matvalues, references, SIMPLIFY = FALSE)
	return(retval)
}

generateSimpleMatrixList <- function(mxModel) {
	retval <- generateMatrixList(mxModel)
	retval <- lapply(retval, function(x) { 
		c(list(as.matrix(x[[1]])), x[-1]) }) 
	return(retval)
}

generateAlgebraList <- function(mxModel) {
	mNames <- names(mxModel@matrices)
	aNames <- names(mxModel@algebras)
	oNames <- names(mxModel@objectives)
	retval <- lapply(mxModel@algebras, generateAlgebraHelper, 
		mNames, append(aNames, oNames))
    return(retval)
}

generateParameterList <- function(flatModel) {
	result <- list()
	if (length(flatModel@matrices) == 0) {
		return(result)
	}
	for(i in 1:length(flatModel@matrices)) {
		matrix <- flatModel@matrices[[i]]
		result <- generateParameterListHelper(
			matrix, result, i - 1L)
	}
	return(result)
}

generateDefinitionList <- function(flatModel) {
	result <- list()
	defLocations <- generateDefinitionLocations(flatModel@datasets)
	if (length(flatModel@matrices) == 0) {
		return(result)
	}
	for(i in 1:length(flatModel@matrices)) {
		result <- generateDefinitionListHelper(
			flatModel@matrices[[i]], 
			result, defLocations, i - 1L)
	}	
	return(result)
}

generateValueList <- function(mxModel, mList, pList) {
	mList <- lapply(mList, function(x) { x[[1]] })
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
	row <- triple[2] + 1
	col <- triple[3] + 1
	return(mList[[mat]][row,col])
}

getObjectiveIndex <- function(flatModel) {
	objective <- flatModel@objective
	if(is.null(objective)) {
		return(NULL)
	} else {
		return(omxLocateIndex(flatModel, objective@name, flatModel@name))
	}
}

updateModelValues <- function(model, flatModel, pList, values) {
	if(length(pList) != length(values)) {
		stop(paste("This model has", length(pList), 
			"parameters, but you have given me", length(values),
			"values"))
	}
	if (length(pList) == 0) {
		return(model)
	}
	for(i in 1:length(pList)) {
		parameters <- pList[[i]]
		parameters <- parameters[3:length(parameters)] # Remove (min, max) bounds
		model <- updateModelValuesHelper(parameters, values[[i]], flatModel, model)
    }
	return(model)
}

updateModelValuesHelper <- function(triples, value, flatModel, model) {
	for(i in 1:length(triples)) {
		triple <- triples[[i]]
		mat <- triple[1] + 1
		row <- triple[2] + 1
		col <- triple[3] + 1
		name <- flatModel@matrices[[mat]]@name
		model[[name]]@values[row,col] <- value
	}
	return(model)
}

removeTail <- function(lst, tailSize) {
    newEnd <- length(lst) - tailSize
    if (newEnd == 0) {
        return(list())
    } else {
        return(lst[1 : newEnd])
    }
}

updateModelMatrices <- function(model, flatModel, values) {
    flatModel@matrices <- removeTail(flatModel@matrices, length(flatModel@constMatrices))
    flatModel@matrices <- removeTail(flatModel@matrices, length(flatModel@freeMatrices))
    values <- removeTail(values, length(flatModel@constMatrices))    
    values <- removeTail(values, length(flatModel@freeMatrices))
	mList <- names(flatModel@matrices)
	if (length(mList) != length(values)) {
		stop(paste("This model has", length(mList), 
			"matrices, but you have given me", length(values),
			"values"))
	}
	if (length(mList) == 0) {
		return(model)
	}	
	model <- updateModelMatricesHelper(mList, values, model)
	return(model)
}


updateModelAlgebras <- function(model, flatModel, values) {
	aNames <- names(flatModel@algebras)
	oNames <- names(flatModel@objectives)
	aList <- append(aNames, oNames)
	if(length(aList) != length(values)) {
		stop(paste("This model has", length(aList), 
			"algebras, but you have given me", length(values),
			"values"))
	}
	if (length(aList) == 0) {
		return(model)
	}	
	model <- updateModelAlgebrasHelper(aList, values, model)
	return(model)
}

updateModelMatricesHelper <- function(mList, values, model) {
	for(i in 1:length(mList)) {
		name <- mList[[i]]
		model[[name]]@values <- values[[i]]
	}
	return(model)
}

updateModelAlgebrasHelper <- function(aList, values, model) {
	for(i in 1:length(aList)) {
		name <- aList[[i]]
		candidate <- model[[name]]
		if (!is.null(candidate) && (length(values[[i]]) > 0) 
			&& !is.nan(values[[i]]) && 
			(is(candidate,"MxAlgebra") || (is(candidate,"MxObjective")))) {
			model[[name]]@result <- as.matrix(values[[i]])
			if (is(candidate, "MxAlgebra")) {
				dimnames(model[[name]]@result) <- dimnames(model[[name]])
			}
		}
	}
	return(model)
}

omxLocateIndex <- function(model, name, referant) {
	if (is.na(name)) { return(as.integer(name)) }
	mNames <- names(model@matrices)
	aNames <- names(model@algebras)
	oNames <- names(model@objectives)
	dNames <- names(model@datasets)		
	matrixNumber <- match(name, mNames)
	algebraNumber <- match(name, append(aNames, oNames))
	dataNumber <- match(name, dNames)
	if (is.na(matrixNumber) && is.na(algebraNumber) && is.na(dataNumber)) {
		msg <- paste("The reference", omxQuotes(name),
			"does not exist.  It is used by the named entity",
			omxQuotes(referant),".")
		stop(msg, call.=FALSE)
	} else if (!is.na(matrixNumber)) {
		return(- matrixNumber)
	} else if (!is.na(dataNumber)) {
		return(dataNumber - 1L)
	} else {
		return(algebraNumber - 1L)
	}
}

omxCheckMatrices <- function(model) {
	matrices <- model@matrices
	lapply(matrices, omxVerifyMatrix)
	if (length(model@submodels) > 0) {
		for(i in 1:length(model@submodels)) {
			submodel <- model@submodels[[i]]
			if(submodel@independent == FALSE) {
				omxCheckMatrices(submodel)
			}
		}
	}
}
