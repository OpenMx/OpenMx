#
#   Copyright 2007-2012 The OpenMx Project
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


generateMatrixList <- function(model) {
	matvalues <- lapply(model@matrices, generateMatrixValuesHelper)
	matnames  <- names(model@matrices)
	names(matvalues) <- matnames
	references <- generateMatrixReferences(model)
	retval <- mapply(function(x,y) { c(list(x), y) }, 
			matvalues, references, SIMPLIFY = FALSE)
	return(retval)
}

generateAlgebraList <- function(model) {
	mNames <- names(model@matrices)
	aNames <- append(names(model@algebras), names(model@objectives))	
	mNumbers <- as.list(as.integer(-1 : (-length(mNames))))
	aNumbers <- as.list(as.integer(0 : (length(aNames) - 1)))
	names(mNumbers) <- mNames
	names(aNumbers) <- aNames
	retval <- lapply(model@algebras, generateAlgebraHelper, 
		mNumbers, aNumbers)
    return(retval)
}

findDependencies <- function(triple, flatModel, dependencies) {
	mNames <- names(flatModel@matrices)
	matrixNum <- triple[[1]] + 1
	matrixName <- mNames[[matrixNum]]
	return(dependencies[[matrixName]])	
}

parameterDependencyList <- function(pList, flatModel, dependencies) {
	if (length(pList) == 2) {
		retval <- list(pList[[1]], pList[[2]], integer())
		return(retval)
	}
	locations <- pList[3:length(pList)]
	deps <- lapply(locations, findDependencies, flatModel, dependencies)
	depnames <- Reduce(union, deps, character())
	depnumbers <- sapply(depnames, doLocateIndex, flatModel, flatModel@name, USE.NAMES=FALSE)
	depnumbers <- as.integer(depnumbers)
	retval <- list(pList[[1]], pList[[2]], depnumbers)
	return(append(retval, locations))
}

generateParameterList <- function(flatModel, dependencies) {
	result <- list()
	if (length(flatModel@matrices) == 0) {
		return(result)
	}
	for(i in 1:length(flatModel@matrices)) {
		matrix <- flatModel@matrices[[i]]
		result <- generateParameterListHelper(
			matrix, result, i - 1L)
	}
	result <- lapply(result, parameterDependencyList, flatModel, dependencies)
	return(result)
}

generateDefinitionList <- function(flatModel, dependencies) {
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
	result <- lapply(result, parameterDependencyList, flatModel, dependencies)
	return(result)
}

generateValueList <- function(mList, pList) {
	mList <- lapply(mList, function(x) { x[[1]] })
	retval <- vector()
	if (length(pList) == 0) {
		return(retval)
	}
	for(i in 1:length(pList)) {
		parameter <- pList[[i]]
		parameter <- parameter[4:length(parameter)] # Remove (min, max, dependencies)
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
		return(imxLocateIndex(flatModel, objective@name, flatModel@name))
	}
}

imxUpdateModelValues <- function(model, flatModel, pList, values) {
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
		parameters <- parameters[4:length(parameters)] # Remove min, max, and dependencies
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
	mList <- names(flatModel@matrices)
	if (length(mList) != length(values)) {
		stop(paste("This model has", length(mList), 
			"matrices, but the backend has returned", length(values),
			"values"))
	}
	if (length(values) == 0) {
		return(model)
	}
	model <- updateModelEntitiesHelper(mList, values, model)
	return(model)
}


updateModelAlgebras <- function(model, flatModel, values) {
	aNames <- names(flatModel@algebras)
	oNames <- names(flatModel@objectives)
	aList <- append(aNames, oNames)
	if(length(aList) != length(values)) {
		stop(paste("This model has", length(aList), 
			"algebras, but the backend has returned", length(values),
			"values"))
	}
	if (length(aList) == 0) {
		return(model)
	}
	model <- updateModelEntitiesHelper(aList, values, model)
	return(model)
}


updateModelEntitiesHelper <- function(entNames, values, model) {
	modelNameMapping <- sapply(entNames, getModelNameString)
	modelNames <- unique(modelNameMapping)
	for(j in 1:length(modelNames)) {
		nextName <- modelNames[[j]]
		selectEnt <- entNames[modelNameMapping == nextName]
		selectVal <- values[modelNameMapping == nextName]
		submodel <- model[[nextName]]
		for(i in 1:length(selectEnt)) {
			name <- selectEnt[[i]]
			candidate <- submodel[[name]]
			value <- selectVal[[i]]
			if (!is.null(candidate) && (length(value) > 0)
				&& !is.nan(value)) {
				if (is(candidate,"MxAlgebra") || is(candidate,"MxObjective")) {
					if (is(candidate, "MxAlgebra")) {
						dimnames(value) <- dimnames(candidate)
						candidate@result <- value
					} else {
						candidate <- objectiveReadAttributes(candidate, value)
					}
				} else if(is(candidate, "MxMatrix")) {
					dimnames(value) <- dimnames(candidate)
					candidate@values <- value
				}
			}
			submodel[[name]] <- candidate
		}
		model[[nextName]] <- submodel
	}
	return(model)
}

imxLocateIndex <- function(model, name, referant) {
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

imxPreprocessModel <- function(model) {
	model@matrices <- lapply(model@matrices, findSquareBrackets)
	model@submodels <- lapply(model@submodels, imxPreprocessModel)
	return(model)
}


imxCheckMatrices <- function(model) {
	matrices <- model@matrices
	lapply(matrices, imxVerifyMatrix)
	submodels <- imxDependentModels(model)
	lapply(submodels, imxCheckMatrices)
}
