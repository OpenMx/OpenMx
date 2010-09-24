#
#   Copyright 2007-2010 The OpenMx Project
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


setClass(Class = "MxRAMObjective",
	representation = representation(
		A = "MxCharOrNumber",
		S = "MxCharOrNumber",
		F = "MxCharOrNumber",
		M = "MxCharOrNumber",
		thresholds = "MxCharOrNumber",
		dims = "character",
		vector = "logical",
		depth = "integer"),
	contains = "MxBaseObjective")

setMethod("initialize", "MxRAMObjective",
	function(.Object, A, S, F, M, dims, thresholds, vector,
		data = as.integer(NA), name = 'objective') {
		.Object@name <- name
		.Object@A <- A
		.Object@S <- S
		.Object@F <- F
		.Object@M <- M
		.Object@data <- data
		.Object@dims <- dims
		.Object@thresholds <- thresholds
		.Object@vector <- vector
		return(.Object)
	}
)

setMethod("genericObjDependencies", signature("MxRAMObjective"),
	function(.Object, dependencies) {
	sources <- c(.Object@A, .Object@S, .Object@F, .Object@M, .Object@thresholds)
	sources <- sources[!is.na(sources)]
	dependencies <- omxAddDependency(sources, .Object@name, dependencies)
	return(dependencies)
})

setMethod("genericObjFunNamespace", signature("MxRAMObjective"), 
	function(.Object, modelname, namespace) {
		.Object@name <- omxIdentifier(modelname, .Object@name)
		.Object@A <- omxConvertIdentifier(.Object@A, modelname, namespace)
		.Object@S <- omxConvertIdentifier(.Object@S, modelname, namespace)
		.Object@F <- omxConvertIdentifier(.Object@F, modelname, namespace)
		.Object@M <- omxConvertIdentifier(.Object@M, modelname, namespace)
		.Object@data <- omxConvertIdentifier(.Object@data, modelname, namespace)
		.Object@thresholds <- sapply(.Object@thresholds, 
			omxConvertIdentifier, modelname, namespace)
		return(.Object)
})

setMethod("genericObjRename", signature("MxRAMObjective"),
	function(.Object, oldname, newname) {
		.Object@A <- renameReference(.Object@A, oldname, newname)
		.Object@S <- renameReference(.Object@S, oldname, newname)
		.Object@F <- renameReference(.Object@F, oldname, newname)
		.Object@M <- renameReference(.Object@M, oldname, newname)
		.Object@data <- renameReference(.Object@data, oldname, newname)
		.Object@thresholds <- sapply(.Object@thresholds, renameReference, oldname, newname)		
		return(.Object)
})

setMethod("genericObjFunConvert", signature("MxRAMObjective", "MxFlatModel"), 
	function(.Object, flatModel, model, defVars) {
		modelname <- omxReverseIdentifier(model, .Object@name)[[1]]	
		name <- .Object@name
		aMatrix <- .Object@A
		sMatrix <- .Object@S
		fMatrix <- .Object@F
		mMatrix <- .Object@M
		data <- .Object@data
		if(is.na(data)) {
			msg <- paste("The RAM objective",
				"does not have a dataset associated with it in model",
				omxQuotes(modelname))
			stop(msg, call. = FALSE)
		}
		mxDataObject <- flatModel@datasets[[.Object@data]]
		if(!is.na(mMatrix) && single.na(mxDataObject@means)) {
			msg <- paste("The RAM objective",
				"has an expected means vector but",
				"no observed means vector in model",
				omxQuotes(modelname))
			stop(msg, call. = FALSE)
		}
		if(!single.na(mxDataObject@means) && is.null(flatModel[[mMatrix]])) {
			msg <- paste("The RAM objective",
				"has an observed means vector but",
				"no expected means vector in model",
				omxQuotes(modelname))
			stop(msg, call. = FALSE)		
		}
		checkNumericData(mxDataObject)
		.Object@A <- omxLocateIndex(flatModel, aMatrix, name)
		.Object@S <- omxLocateIndex(flatModel, sMatrix, name)
		.Object@F <- omxLocateIndex(flatModel, fMatrix, name)
		.Object@M <- omxLocateIndex(flatModel, mMatrix, name)
		.Object@data <- as.integer(omxLocateIndex(flatModel, data, name))
		verifyObservedNames(mxDataObject@observed, mxDataObject@type, flatModel, modelname, "RAM")
		fMatrix <- flatModel[[fMatrix]]@values
		if (is.null(dimnames(fMatrix))) {
			msg <- paste("The F matrix of model",
				omxQuotes(modelname), "does not contain dimnames")
			stop(msg, call. = FALSE)		
		}
		if (is.null(dimnames(fMatrix)[[2]])) {
			msg <- paste("The F matrix of model",
				omxQuotes(modelname), "does not contain colnames")
			stop(msg, call. = FALSE)
		}
		mMatrix <- flatModel[[mMatrix]]		
		if(!is.null(mMatrix)) {
			means <- dimnames(mMatrix@values)
			if (is.null(means)) {
				msg <- paste("The M matrix associated",
				"with the RAM objective function in model", 
				omxQuotes(modelname), "does not contain dimnames.")
				stop(msg, call.=FALSE)	
			}
			meanRows <- means[[1]]
			meanCols <- means[[2]]
			if (!is.null(meanRows) && length(meanRows) > 1) {
				msg <- paste("The M matrix associated",
				"with the RAM objective in model", 
				omxQuotes(modelname), "is not a 1 x N matrix.")
				stop(msg, call.=FALSE)
			}
			if (!identical(dimnames(fMatrix)[[2]], meanCols)) {
				msg <- paste("The column names of the F matrix",
					"and the column names of the M matrix",
					"in model", 
					omxQuotes(modelname), "do not contain identical",
					"names.")
				stop(msg, call.=FALSE)
			}
		}
		translatedNames <- fMatrixTranslateNames(fMatrix, modelname)
		if (!identical(translatedNames, rownames(mxDataObject@observed))) {
			msg <- paste("The names of the manifest",
				"variables in the F matrix of model",
				omxQuotes(modelname), "does not match the",
				"dimnames of the observed covariance matrix")
			stop(msg, call. = FALSE)
		}
		.Object@depth <- generateRAMDepth(flatModel, aMatrix, model@options)
		return(.Object)
})

generateRAMDepth <- function(flatModel, aMatrixName, modeloptions) {
	mxObject <- flatModel[[aMatrixName]]
	if (!is(mxObject, "MxMatrix")) {
		return(as.integer(NA))
	}
	if (identical(modeloptions[['RAM Inverse Optimization']], "No")) {
		return(as.integer(NA))
	}
	if (is.null(modeloptions[['RAM Inverse Optimization']]) &&
		identical(getOption('mxOptions')[['RAM Inverse Optimization']], "No")) {
		return(as.integer(NA))
	}	
	maxdepth <- modeloptions[['RAM Max Depth']]
	if (is.null(maxdepth) || (length(maxdepth) != 1) ||
		is.na(maxdepth) || !is.numeric(maxdepth) || maxdepth < 0) {
		maxdepth <- nrow(mxObject) - 1
	}
	return(omxGetRAMDepth(mxObject, maxdepth))
}

omxGetRAMDepth <- function(A, maxdepth = nrow(A) - 1) {
	mxObject <- A
	aValues <- matrix(0, nrow(mxObject), ncol(mxObject))
	aValues[mxObject@free] <- 1
	aValues[mxObject@values != 0] <- 1
	return(generateDepthHelper(aValues, aValues, 0, maxdepth))
}

generateDepthHelper <- function(aValues, currentProduct, depth, maxdepth) {
	if (depth > maxdepth) {
		return(as.integer(NA))
	}
	if (all(currentProduct == 0)) { 
		return(as.integer(depth))
	} else {
		return(generateDepthHelper(aValues, currentProduct %*% aValues, depth + 1, maxdepth))
	}
}

fMatrixTranslateNames <- function(fMatrix, modelName) {
	retval <- character()
	colNames <- dimnames(fMatrix)[[2]]
	for(i in 1:nrow(fMatrix)) {
		irow <- fMatrix[i,]
		matches <- which(irow == 1)
		if (length(matches) != 1) {
			err <- paste("The model",
				omxQuotes(modelName), "does not contain",
				"a valid F matrix")
			stop(err, call. = FALSE)
		}
		retval[[i]] <- colNames[[matches[[1]]]]
	}
	return(retval)
}

updateRAMdimnames <- function(flatObjective, job, modelname) {
	fMatrixName <- flatObjective@F
	mMatrixName <- flatObjective@M
	if (is.na(mMatrixName)) {
		mMatrix <- NA
	} else {
		mMatrix <- job[[mMatrixName]]
	}
	fMatrix <- job[[fMatrixName]]
	if (is.null(fMatrix)) {
		stop(paste("Unknown F matrix name", 
			omxQuotes(simplifyName(fMatrixName, modelname)),
			"detected in the objective function",
			"of model", omxQuotes(modelname)), call. = FALSE)
	}
	dims <- flatObjective@dims
	if (!is.null(dimnames(fMatrix)) && !single.na(dims) && 
		!identical(dimnames(fMatrix)[[2]], dims)) {
		msg <- paste("The F matrix associated",
			"with the RAM objective in model", 
			omxQuotes(modelname), "contains dimnames and",
			"the objective function has specified dimnames")
		stop(msg, call.=FALSE)		
	}
	if (is.null(dimnames(fMatrix)) && !single.na(dims)) {
		dimnames(fMatrix) <- list(c(), dims)
		job[[fMatrixName]] <- fMatrix
	}
	if (!isS4(mMatrix) && (is.null(mMatrix) || is.na(mMatrix))) return(job)
	if (!is.null(dimnames(mMatrix)) && !single.na(dims) &&
		!identical(dimnames(mMatrix), list(NULL, dims))) {
		msg <- paste("The M matrix associated",
			"with the RAM objective in model", 
			omxQuotes(modelname), "contains dimnames and",
			"the objective function has specified dimnames")
		stop(msg, call.=FALSE)	
	}
	if (is.null(dimnames(mMatrix)) && !single.na(dims)) {
		dimnames(mMatrix) <- list(NULL, dims)
		job[[mMatrixName]] <- mMatrix
	}
	return(job)
}

setMethod("genericObjModelConvert", "MxRAMObjective",
	function(.Object, job, model, namespace, flatJob) {
		if(is.na(.Object@data)) {
			msg <- paste("The RAM objective",
				"does not have a dataset associated with it in model",
				omxQuotes(model@name))
			stop(msg, call.=FALSE)
		}
		job <- updateRAMdimnames(.Object, job, model@name)
		if (flatJob@datasets[[.Object@data]]@type != 'raw') {
			if (.Object@vector) {
				msg <- paste("The RAM objective",
					"in model", omxQuotes(model@name), "has specified",
					"'vector' = TRUE, but the observed data is not raw data")
				stop(msg, call.=FALSE)
			}
			return(job)
		}
		if (is.na(.Object@M) || is.null(flatJob[[.Object@M]])) {
			msg <- paste("The RAM objective",
				"has raw data but is missing",
				"an expected means vector in model",
				omxQuotes(model@name))
			stop(msg, call.=FALSE)			
		}
		dims <- dimnames(flatJob[[.Object@F]])
		if (is.null(dims) || is.null(dims[[2]])) {
			msg <- paste("The F matrix in model",
				omxQuotes(model@name),
				"does not have column names.")
			stop(msg, call.=FALSE)
		}
		if (availableName(model, namespace, 'I')) {
			iName <- 'I'
		} else {
			iName <- omxUntitledName()
		}
		model <- mxModel(model, mxMatrix(type="Iden", 
			nrow=nrow(flatJob[[.Object@A]]), name = iName))
		if (availableName(model, namespace, 'Z')) {
			zName <- 'Z'
		} else {
			zName <- omxUntitledName()
		}
		zFormula <- substitute(solve(I - A),
			list(I = as.symbol(iName), A = as.symbol(.Object@A)))
		algebra <- eval(substitute(mxAlgebra(x, y),
			list(x = zFormula, y = zName)))
		model <- mxModel(model, algebra)
		if (availableName(model, namespace, 'covariance')) {
			covName <- 'covariance'
		} else {
			covName <- omxUntitledName()
		}
		covFormula <- substitute(F %*% Z %*% S %*% t(Z) %*% t(F),
			list(F = as.symbol(.Object@F), Z = as.symbol(zName),
				S = as.symbol(.Object@S)))
		algebra <- eval(substitute(mxAlgebra(x, y),
			list(x = covFormula, y = covName)))
		translatedNames <- fMatrixTranslateNames(flatJob[[.Object@F]]@values, model@name)
		dimnames(algebra) <- list(translatedNames, translatedNames)
		model <- mxModel(model, algebra)
		meansFormula <- substitute(t(F %*% Z %*% t(M)),
			list(F = as.symbol(.Object@F), Z = as.symbol(zName),
				M = as.symbol(.Object@M)))
		if (availableName(model, namespace, 'means')) {
			meansName <- 'means'
		} else {
			meansName <- omxUntitledName()
		}
		algebra <- eval(substitute(mxAlgebra(x, y),
			list(x = meansFormula, y = meansName)))
		dimnames(algebra) <- list(NULL, translatedNames)
		model <- mxModel(model, algebra)
		objective <- eval(substitute(mxFIMLObjective(covariance = x, 
			means = y, thresholds = z, vector = w),
			list(x = covName, y = meansName, z = .Object@thresholds, w = .Object@vector)))
		objective@.translated <- TRUE
		metadata <- new("MxRAMMetaData", .Object@A, .Object@S, .Object@F, 
			.Object@M, generateRAMDepth(flatJob, .Object@A, job@options))
		objective@metadata <- metadata
		model@objective <- objective
		class(model) <- 'MxModel'
		job[[model@name]] <- model
		return(job)
	}
)

mxRAMObjective <- function(A, S, F, M = NA, dimnames = NA, thresholds = NA, vector = FALSE) {
	if (missing(A) || typeof(A) != "character") {
		msg <- paste("argument 'A' is not a string",
			"(the name of the 'A' matrix)")
		stop(msg)
	}	
	if (missing(S) || typeof(S) != "character") {
		msg <- paste("argument 'S' is not a string",
			"(the name of the 'S' matrix)")
		stop(msg)
	}
	if (missing(F) || typeof(F) != "character") {
		msg <- paste("argument 'F' is not a string",
			"(the name of the 'F' matrix)")
		stop(msg)
	}
	if (!(single.na(M) || typeof(M) == "character")) {
		msg <- paste("argument M is not a string",
			"(the name of the 'M' matrix)")
		stop(msg)
	}
	if (is.na(M)) M <- as.integer(NA)
	if (single.na(thresholds)) thresholds <- as.character(NA)
	if (single.na(dimnames)) dimnames <- as.character(NA)
	if (!is.vector(dimnames) || typeof(dimnames) != 'character') {
		stop("Dimnames argument is not a character vector")
	}
	if (length(thresholds) != 1) {
		stop("Thresholds argument must be a single matrix or algebra name")
	}
	if (length(dimnames) == 0) {
		stop("Dimnames argument cannot be an empty vector")
	}
	if (length(dimnames) > 1 && any(is.na(dimnames))) {
		stop("NA values are not allowed for dimnames vector")
	}
	if (length(vector) > 1 || typeof(vector) != "logical") {
		stop("Vector argument is not a logical value")
	}
	return(new("MxRAMObjective", A, S, F, M, dimnames, thresholds, vector))
}

displayRAMObjective <- function(objective) {
	cat("MxRAMObjective", omxQuotes(objective@name), '\n')
	cat("@A :", omxQuotes(objective@A), '\n')
	cat("@S :", omxQuotes(objective@S), '\n')
	cat("@F :", omxQuotes(objective@F), '\n')
	if (is.na(objective@M)) {
		cat("@M :", objective@M, '\n')
	} else {
		cat("@M :", omxQuotes(objective@M), '\n')
	}
	cat("@vector :", objective@vector, '\n')
	if (single.na(objective@dims)) {
		cat("@dims : NA \n")
	} else {
		cat("@dims :", omxQuotes(objective@dims), '\n')
	}		
	if (single.na(objective@thresholds)) {
		cat("@thresholds : NA \n")
	} else {
		cat("@thresholds :", omxQuotes(objective@thresholds), '\n')
	}
	if (length(objective@result) == 0) {
		cat("@result: (not yet computed) ")
	} else {
		cat("@result:\n")
	}
	print(objective@result)
	invisible(objective)
}

setMethod("print", "MxRAMObjective", function(x,...) { 
	displayRAMObjective(x) 
})

setMethod("show", "MxRAMObjective", function(object) { 
	displayRAMObjective(object) 
})
