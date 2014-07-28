#
#   Copyright 2007-2014 The OpenMx Project
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
	aNames <- append(names(model@algebras), names(model@fitfunctions))	
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

	
isExpectation <- function(name) {
	return(length(grep("expectation", name, fixed=TRUE)) > 0)
}

parameterDependencyList <- function(pList, flatModel, dependencies) {
	if (length(pList) == 3) {
		return(retval)
	}
	locations <- pList[4:length(pList)]
	deps <- lapply(locations, findDependencies, flatModel, dependencies)
	depnames <- Reduce(union, deps, character())
	depnames <- Filter(Negate(isExpectation), depnames)
	depnumbers <- sapply(depnames, doLocateIndex, flatModel, flatModel@name, USE.NAMES=FALSE)
	depnumbers <- as.integer(depnumbers)
	retval <- append(pList[1:3], list(depnumbers))
	append(retval, locations)
}

buildFreeVarGroupList <- function(flatModel) {
	if (is.null(flatModel@compute)) return()
	got <- getFreeVarGroup(flatModel@compute)
	if (length(got)) {
		if (anyDuplicated(unlist(got[seq(1,length(got),2)]))) {
			stop("freeSet IDs are not unique") #for debugging
		}
		members <- unlist(got[seq(2,length(got),2)])
		recog <- match(members, names(flatModel@matrices))
		if (any(is.na(recog))) {
			stop(paste("freeSet component", omxQuotes(members[is.na(recog)]), "not found"))
		}
	}
	got
}

generateParameterList <- function(flatModel, dependencies, freeVarGroups) {
	mList <- flatModel@matrices
	pList <- list()
	for(i in 1:length(mList)) {
		matrix <- mList[[i]]
		pList <- generateParameterListHelper(matrix, pList, i - 1L, freeVarGroups)
	}
	pList <- lapply(pList, parameterDependencyList, flatModel, dependencies)

	if (length(pList)) for(i in 1:length(pList)) {
		original <- pList[[i]]
		svalues <- original[5:length(original)]
		svalue <- NA
		if (length(svalues) > 1) {
			values <- sapply(svalues, generateValueHelper, mList)
			if (!all(values == values[[1]])) {
				warning(paste('Parameter',i,'has multiple start values.',
					      'Selecting', values[[1]]))
			}
			svalue <- values[[1]]
		} else {
			svalue <- generateValueHelper(svalues[[1]], mList)
		}
		pList[[i]] <- c(original, svalue)
	}
	flatModel@parameters <- pList
	flatModel
}

definitionDependencyList <- function(pList, flatModel, dependencies) {
	if (length(pList) == 2) {
		retval <- list(pList[[1]], pList[[2]], integer())
		return(retval)
	}
	locations <- pList[3:length(pList)]
	deps <- lapply(locations, findDependencies, flatModel, dependencies)
	depnames <- Reduce(union, deps, character())
	depnames <- Filter(Negate(isExpectation), depnames)
	depnumbers <- sapply(depnames, doLocateIndex, flatModel, flatModel@name, USE.NAMES=FALSE)
	depnumbers <- as.integer(depnumbers)
	retval <- list(pList[[1]], pList[[2]], depnumbers)
	return(append(retval, locations))
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
	result <- lapply(result, definitionDependencyList, flatModel, dependencies)
	return(result)
}

generateValueHelper <- function(triple, mList) {
	mat <- triple[1] + 1
	row <- triple[2] + 1
	col <- triple[3] + 1
	val <- mList[[mat]]@values[row,col]
	if (is.na(val)) {
		stop(paste("Starting value in ",names(mList)[[mat]],
			   "[",row,",",col,"] is missing", sep=""))
	}
	return(val)
}

imxUpdateModelValues <- function(model, flatModel, values) {
	pList <- flatModel@parameters
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
		parameters <- parameters[5:(length(parameters)-1)]
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
	oNames <- names(flatModel@fitfunctions)
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

updateModelExpectations <- function(model, flatModel, values) {
	eNames <- names(flatModel@expectations)
	if (length(eNames) != length(values)) {
		stop(paste("This model has", length(eNames),
			   "expectations, but the backend has returned", length(values),
			   "values"))
	}
	if (length(eNames) == 0) return(model)
	updateModelEntitiesHelper(eNames, values, model)
}

updateModelData <- function(model, flatModel, values) {
	dNames <- names(flatModel@datasets)
	if (length(dNames) != length(values)) {
		stop(paste("This model has", length(dNames),
			   "data, but the backend has returned", length(values),
			   "values"))
	}
	if (length(dNames) == 0) return(model)
	updateModelEntitiesHelper(dNames, values, model)
}

updateModelEntitiesTargetModel <- function(model, entNames, values, modelNameMapping) {
    nextName <- model@name
    selectEnt <- entNames[modelNameMapping == nextName]
    selectVal <- values[modelNameMapping == nextName]
    if (length(selectEnt) > 0) {
		for(i in 1:length(selectEnt)) {
			name <- selectEnt[[i]]
			candidate <- model[[name]]
			value <- selectVal[[i]]
			if (!is.null(candidate) && (length(value) > 0)
				&& !is.nan(value)) {
				if (is(candidate, "MxAlgebra")) {
					cdim <- sapply(dimnames(candidate), length)
					mask <- cdim != 0
					if (any(cdim[mask] != dim(value)[mask])) {
						warning(paste(paste(model@name, candidate@name, sep="."),
							      "is dimension", omxQuotes(dim(value)),
							      "but the dimnames imply dimension",
							      omxQuotes(cdim)))
					} else {
						dimnames(value) <- dimnames(candidate)
					}
					candidate@result <- value
				} else if (is(candidate,"MxFitFunction")) {
					candidate@result <- as.matrix(value)
					attr <- attributes(value)
					attr$dim <- NULL
					candidate@info <- attr
				} else if(is(candidate, "MxMatrix")) {
					if (any(dim(candidate@values) != dim(value))) {
						msg <- paste("Backend returned a", omxQuotes(dim(value)),
							     "matrix for a", omxQuotes(dim(candidate@values)),
							     "matrix. Not sure how to proceed.",
							     "Please report this error to the OpenMx support team.")
						stop(msg)
					} else {
						dimnames(value) <- dimnames(candidate)
						candidate@values <- value
					}
				} else if (is(candidate, "MxExpectation")) {
					for (sl in names(attributes(value))) {
						slot(candidate, sl) <- attr(value, sl)
					}
				} else if (is(candidate, "MxDataDynamic")) {
					candidate@numObs <- value
				}
				model[[name]] <- candidate
			}
		}
	}
    if (length(model@submodels) > 0) {
        model@submodels <- lapply(model@submodels, 
            updateModelEntitiesTargetModel, entNames, values, modelNameMapping)
    }
	return(model)
}

updateModelEntitiesHelper <- function(entNames, values, model) {
    modelNameMapping <- sapply(entNames, getModelNameString)
    model <- updateModelEntitiesTargetModel(model, entNames, 
		values, modelNameMapping)
    return(model)
}

clearModifiedSinceRunRecursive <- function(model) {
	if (length(model@submodels) > 0) {
		model@submodels <- lapply(model@submodels, clearModifiedSinceRunRecursive)
	}
	model@.modifiedSinceRun <- FALSE
	model
}

##' imxLocateIndex
##'
##' This is an internal function exported for those people who know
##' what they are doing.
##'
##' @param model model
##' @param name name
##' @param referant referant
imxLocateIndex <- function(model, name, referant) {
	if (length(name) == 0) return(name)
	if (is.na(name)) { return(as.integer(name)) }
	mNames <- names(model@matrices)
	aNames <- names(model@algebras)
	fNames <- names(model@fitfunctions)
	eNames <- names(model@expectations)
	dNames <- names(model@datasets)		
	matrixNumber <- match(name, mNames)
	algebraNumber <- match(name, append(aNames, fNames))
	dataNumber <- match(name, dNames)
	expectationNumber <- match(name, eNames)
	if (is.na(matrixNumber) && is.na(algebraNumber) 
		&& is.na(dataNumber) && is.na(expectationNumber)) {
		reftype <- "named reference"
		if (typeof(referant) == "S4") {
			reftype <- 'object'
		}
		msg <- paste("The reference", omxQuotes(name),
			"does not exist.  It is used by",
			     reftype, omxQuotes(referant), ".")
		stop(msg, call.=FALSE)
	} else if (!is.na(matrixNumber)) {
		return(- matrixNumber)
	} else if (!is.na(dataNumber)) {
		return(dataNumber - 1L)
	} else if (!is.na(expectationNumber)) {
		return(expectationNumber - 1L)
	} else {
		return(algebraNumber - 1L)
	}
}

##' imxPreprocessModel
##'
##' This is an internal function exported for those people who know
##' what they are doing.
##'
##' @param model model
imxPreprocessModel <- function(model) {
	model@matrices <- lapply(model@matrices, findSquareBrackets)
	model@submodels <- lapply(model@submodels, imxPreprocessModel)
	return(model)
}


##' imxCheckMatrices
##'
##' This is an internal function exported for those people who know
##' what they are doing.
##'
##' @param model model
imxCheckMatrices <- function(model) {
	matrices <- model@matrices
	lapply(matrices, imxVerifyMatrix)
	submodels <- imxDependentModels(model)
	lapply(submodels, imxCheckMatrices)
}
