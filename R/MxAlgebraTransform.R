#
#   Copyright 2007-2016 The OpenMx Project
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

transformAlgebras <- function(model) {
	model <- transformAlgebraHelper(model, model)
}

transformAlgebraHelper <- function(model, job) {
	repeat {
		topBindExpressions <- locateTopBindExpressions(model, job)
		subBindExpressions <- locateSubBindExpressions(model, job)
		if (length(topBindExpressions) + length(subBindExpressions) == 0) break
		topBindMatrices <- generateTopBindMatrices(topBindExpressions, model, job)
		subBindMatrices <- generateSubBindMatrices(subBindExpressions, model, job)
		model <- addBindMatrices(model, topBindMatrices, subBindMatrices)
		model <- replaceBindExpressions(model, subBindExpressions, subBindMatrices)
		if (job@name == model@name) {
			job <- model
		} else {
			job[[model@name]] <- model
		}
	}
	if (length(model@submodels) > 0) {
		for(i in 1:length(model@submodels)) {
			job <- transformAlgebraHelper(model@submodels[[i]], job)
		}
	}
	return(job)
}

addBindMatrices <- function(model, topBindMatrices, subBindMatrices) {
	model@algebras[names(topBindMatrices)] <- NULL
	names(subBindMatrices) <- imxExtractNames(subBindMatrices)
	model@matrices <- c(model@matrices, topBindMatrices, subBindMatrices)
	return(model)
}

replaceBindExpressions <- function(model, bindExpressions, bindMatrices) {
	model@algebras <- lapply(model@algebras, replaceBindHelper, bindExpressions, bindMatrices)
	return(model)
}

replaceBindHelper <- function(algebra, bindExpressions, bindMatrices) {
	algebra@formula <- replaceBindFormula(algebra@formula, bindExpressions, bindMatrices)
	return(algebra)
}

replaceBindFormula <- function(formula, bindExpressions, bindMatrices) {
	if (length(formula) == 1) {
		return(formula)
	}
	if (startsWithBind(formula)) {
		index <- sapply(bindExpressions, identical, formula)
		index <- match(TRUE, index)
	} else {
		index <- NA
	}
	if (!is.na(index)) {
		return(as.symbol(bindMatrices[[index]]@name))
	} else {
		for(i in 2:length(formula)) {
			formula[[i]] <- replaceBindFormula(formula[[i]], bindExpressions, bindMatrices)
		}
		return(formula)
	}
}

locateTopBindExpressions <- function(model, job) {
	retval <- lapply(model@algebras, locateBindHelper, model, job)
	if (length(retval) > 0) {
		retval <- unlist(retval, recursive = TRUE)
	}
	return(retval)
}

locateSubBindExpressions <- function(model, job) {
	retval <- lapply(model@algebras, locateSubBindHelper, model, job)
	if (length(retval) > 0) {
		retval <- unlist(retval, recursive = TRUE)
	}
	return(retval)
}

locateTopBindExpressions <- function(model, job) {
	retval <- list()
	if (length(model@algebras) == 0) return(retval)
	for(i in 1:length(model@algebras)) {
		algebra <- model@algebras[[i]]
		formula <- algebra@formula
		name <- algebra@name
		if (length(formula) > 1 && startsWithBind(formula) && allArgsMatrices(formula, model, job)) {
			retval[[name]] <- formula
		}
	}
	return(retval)
}

locateSubBindHelper <- function(algebra, model, job) {
	return(locateBindFormula(algebra@formula, model, job, TRUE))
}

locateBindFormula <- function(formula, model, job, skip) {
	if (length(formula) == 1) {
		return(list())
	}
	if (startsWithBind(formula) && allArgsMatrices(formula, model, job) && !skip) {
		return(formula)
	} else {
		return(lapply(formula[-1], locateBindFormula, model, job, FALSE))
	}
}

startsWithBind <- function(formula) { 
	return(formula[[1]] == 'cbind' || formula[[1]] == 'rbind')
}

allArgsMatrices <- function(formula, model, job) {
	formula <- formula[-1]
	return(all(sapply(formula, isArgMatrix, model, job)))
}

isArgMatrix <- function(arg, model, job) {
	if (length(arg) > 1) return(FALSE)
	arg <- as.character(arg)
	if (imxIsDefinitionVariable(arg)) return(FALSE)
	components <- unlist(strsplit(arg, imxSeparatorChar, fixed = TRUE))
	if (length(components) == 2) {
		object <- job[[arg]]
		if(is.null(object)) return(FALSE)
		return(is(object, "MxMatrix"))
	} else {
		object <- model[[arg]]
		if(is.null(object)) return(FALSE)
		return(is(object, "MxMatrix"))
	}
}

generateTopBindMatrices <- function(bindExpressions, model, job) {
	retval <- list()
	if (length(bindExpressions) == 0) return(retval)
	for(i in 1:length(bindExpressions)) {
		bindExpression <- bindExpressions[[i]]
		name <- names(bindExpressions)[[i]]
		algebraObject <- model[[name]]
		dimnames <- algebraObject@.dimnames
		matrices <- lapply(bindExpression[-1], findMatrix, model, job)
		free <- lapply(matrices, slot, "free")
		values <- lapply(matrices, slot, "values")
		labels <- lapply(matrices, slot, "labels")
		lbound <- lapply(matrices, slot, "lbound")
		ubound <- lapply(matrices, slot, "ubound")
		functionName <- as.character(bindExpression[[1]])
		free <- do.call(functionName, free)
		values <- do.call(functionName, values)
		labels <- do.call(functionName, labels)
		lbound <- do.call(functionName, lbound)
		ubound <- do.call(functionName, ubound)
		newmatrix <- mxMatrix('Full', nrow(free), ncol(free), free, values,
			labels, lbound, ubound, dimnames = dimnames, name = name)
		retval[[name]] <- newmatrix
	}
	return(retval)
}

generateSubBindMatrices <- function(bindExpressions, model, job) {
	retval <- lapply(bindExpressions, generateSubBindMatrix, model, job)
	return(retval)
}

generateSubBindMatrix <- function(bindExpression, model, job) {
	matrices <- lapply(bindExpression[-1], findMatrix, model, job)
	free <- lapply(matrices, slot, "free")
	values <- lapply(matrices, slot, "values")
	labels <- lapply(matrices, slot, "labels")
	lbound <- lapply(matrices, slot, "lbound")
	ubound <- lapply(matrices, slot, "ubound")
	functionName <- as.character(bindExpression[[1]])
	free <- do.call(functionName, free)
	values <- do.call(functionName, values)
	labels <- do.call(functionName, labels)
	lbound <- do.call(functionName, lbound)
	ubound <- do.call(functionName, ubound)
	retval <- eval(substitute(mxMatrix('Full', free = a, values = b,
		labels = c, lbound = d, ubound = e), list(
		a = free, b = values, c = labels, d = lbound, e = ubound)))
	return(retval)
}

findMatrix <- function(symbol, model, job) {
	name <- as.character(symbol)
	components <- unlist(strsplit(name, imxSeparatorChar, fixed = TRUE))
	if (length(components) == 2) {
		return(job[[name]])
	} else {
		return(model[[name]])
	}
}
