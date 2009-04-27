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

omxCheckNamedEntity <- function(model, slotname, nameList) {
	entity <- slot(model, slotname)
	entityNames <- names(entity)
	omxNameAlignment(entityNames, omxGetNames(entity))
	entityIntersect <- intersect(entityNames, nameList)
	if (length(entityIntersect) > 0) {
		stop(omxNamespaceErrorMessage(entityIntersect), call.=FALSE)
	}
	nameList <- append(nameList, entityNames)
	return(nameList)
}

omxCheckDataColumns <- function(dataset, nameList) {
	columnNames <- dimnames(dataset)[[2]]
	nameIntersect <- intersect(columnNames, nameList)
	if (length(nameIntersect) > 0) {
		stop(omxNamespaceErrorMessage(nameIntersect), call.=FALSE)
	}
	nameList <- append(nameList, columnNames)	
	return(nameList)
}


omxCheckNamespaceHelper <- function(model, nameList) {
	nameList <- omxCheckNamedEntity(model, "matrices", nameList)
	nameList <- omxCheckNamedEntity(model, "algebras", nameList)
	nameList <- omxCheckNamedEntity(model, "submodels", nameList)
	nameList <- omxCheckNamedEntity(model, "constraints", nameList)
	nameList <- omxCheckDataColumns(model@data, nameList)
	if (!is.null(model@objective) && (model@objective@name %in% nameList)) {
		stop(omxNamespaceErrorMessage(model@objective@name), call.=FALSE)
	} else if(!is.null(model@objective)) {
		nameList <- append(nameList, model@objective@name)
	}
	if (!is.null(model@data) && (model@data@name %in% nameList)) {
		stop(omxNamespaceErrorMessage(model@data@name), call.=FALSE)
	} else if(!is.null(model@data)) {
		nameList <- append(nameList, model@data@name)
	}
	if (length(model@submodels) > 0) {
		for(i in 1:length(model@submodels)) {
			nameList <- omxCheckNamespaceHelper(model@submodels[[i]], nameList)
		}
	}
	return(nameList)
}

omxShareData <- function(model) {
	return(omxShareDataHelper(model, model@data))
}

omxShareDataHelper <- function(model, current) {
	if((is.null(model@data)) && (model@independent == TRUE)) {
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
	type <- class(mxMatrix)[[1]]
	rows <- nrow(mxMatrix)
	cols <- ncol(mxMatrix)	
	return(new(type, mxMatrix@name, mxMatrix@values, 
		FALSE, mxMatrix@labels, rows, cols))
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
	defaultData <- model@data
	res <- omxFlattenModelHelper(model, res, defaultData)
	res@submodels <- list()
	return(res)
}

omxFlattenModelHelper <- function(model, dest, defaultData) {
	if (!is.null(model@objective)) {
		if(is.null(model@data) && !is.null(defaultData)) {
			model@objective@data <- defaultData@name			
		} else if (!is.null(model@data)) {
			model@objective@data <- model@data@name
		}
		dest@objectives[[model@objective@name]] <- model@objective
	}
	if (!is.null(model@data)) {
		dest@datasets[[model@data@name]] <- model@data
	}
	if (is.null(defaultData)) {
		defaultData <- model@data
	}
	if (length(model@submodels) > 0) {
		for(i in 1:length(model@submodels)) {
			submodel <- model@submodels[[i]]
			dest@matrices    <- append(dest@matrices, submodel@matrices)
			dest@algebras    <- append(dest@algebras, submodel@algebras)
			dest@constraints <- append(dest@constraints, 
				submodel@constraints) 
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
