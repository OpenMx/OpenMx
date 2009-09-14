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


setClass(Class = "MxRAMObjective",
	representation = representation(
		A = "MxCharOrNumber",
		S = "MxCharOrNumber",
		F = "MxCharOrNumber",
		M = "MxCharOrNumber",
		thresholds = "MxCharOrNumber"),
	contains = "MxBaseObjective")

setMethod("initialize", "MxRAMObjective",
	function(.Object, A, S, F, M, thresholds,  
		data = as.integer(NA), name = 'objective') {
		.Object@name <- name
		.Object@A <- A
		.Object@S <- S
		.Object@F <- F
		.Object@M <- M
		.Object@data <- data
		.Object@thresholds <- thresholds
		.Object@dependencies <- c('A', 'S', 'F', 'M', 'thresholds')
		return(.Object)
	}
)

setMethod("omxObjFunNamespace", signature("MxRAMObjective"), 
	function(.Object, modelname, namespace) {
		.Object@name <- omxIdentifier(modelname, .Object@name)
		.Object@A <- omxConvertIdentifier(.Object@A, modelname, namespace)
		.Object@S <- omxConvertIdentifier(.Object@S, modelname, namespace)
		.Object@F <- omxConvertIdentifier(.Object@F, modelname, namespace)
		.Object@M <- omxConvertIdentifier(.Object@M, modelname, namespace)
		.Object@data <- omxConvertIdentifier(.Object@data, modelname, namespace)
		.Object@thresholds <- omxConvertIdentifier(.Object@thresholds, modelname, namespace)
		return(.Object)
})

setMethod("omxObjFunConvert", signature("MxRAMObjective", "MxFlatModel"), 
	function(.Object, flatModel, model) {
		name <- .Object@name
		aMatrix <- .Object@A
		sMatrix <- .Object@S
		fMatrix <- .Object@F
		mMatrix <- .Object@M
		data <- .Object@data
		if(is.na(data)) {
			msg <- paste("The RAM objective",
				"does not have a dataset associated with it in model",
				omxQuotes(flatModel@name))
			stop(msg, call. = FALSE)
		}
		.Object@A <- omxLocateIndex(flatModel, aMatrix, name)
		.Object@S <- omxLocateIndex(flatModel, sMatrix, name)
		.Object@F <- omxLocateIndex(flatModel, fMatrix, name)
		.Object@M <- omxLocateIndex(flatModel, mMatrix, name)
		.Object@data <- as.integer(omxLocateIndex(flatModel, data, name))
		return(.Object)
})

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

setMethod("omxObjModelConvert", "MxRAMObjective",
	function(.Object, flatModel, model) {
		if(is.na(.Object@data)) {
			msg <- paste("The RAM objective",
				"does not have a dataset associated with it in model",
				omxQuotes(model@name))
			stop(msg, call.=FALSE)
		}		
		if (flatModel@datasets[[.Object@data]]@type != 'raw' || 
			is.na(.Object@M)) {
			return(model)
		}
		dims <- dimnames(flatModel[[.Object@F]])
		if (is.null(dims) || is.null(dims[[2]])) {
			msg <- paste("The F matrix in model",
				omxQuotes(model@name),
				"does not have column names.")
			stop(msg, call.=FALSE)
		}
		if (is.null(model[['I']])) {
			iName <- 'I'
		} else {
			iName <- omxUntitledName()
		}
		model <- mxModel(model, mxMatrix(type="Iden", 
			nrow=nrow(flatModel[[.Object@A]]), name = iName))
		if (is.null(model[['Z']])) {
			zName <- 'Z'
		} else {
			zName <- omxUntitledName()
		}
		zFormula <- substitute(solve(I - A),
			list(I = as.symbol(iName), A = as.symbol(.Object@A)))
		algebra <- eval(substitute(mxAlgebra(x, y),
			list(x = zFormula, y = zName)))
		model <- mxModel(model, algebra)
		if (is.null(model[['covariance']])) {
			covName <- 'covariance'
		} else {
			covName <- omxUntitledName()
		}
		covFormula <- substitute(F %*% Z %*% S %*% t(Z) %*% t(F),
			list(F = as.symbol(.Object@F), Z = as.symbol(zName),
				S = as.symbol(.Object@S)))
		algebra <- eval(substitute(mxAlgebra(x, y),
			list(x = covFormula, y = covName)))
		translatedNames <- fMatrixTranslateNames(flatModel[[.Object@F]]@values, model@name)
		dimnames(algebra) <- list(translatedNames, translatedNames)
		model <- mxModel(model, algebra)
		meansFormula <- substitute(t(F %*% Z %*% t(M)),
			list(F = as.symbol(.Object@F), Z = as.symbol(zName),
				M = as.symbol(.Object@M)))
		if (is.null(model[['means']])) {
			meansName <- 'means'
		} else {
			meansName <- omxUntitledName()
		}
		algebra <- eval(substitute(mxAlgebra(x, y),
			list(x = meansFormula, y = meansName)))
		dimnames(algebra) <- list(NULL, translatedNames)
		model <- mxModel(model, algebra)
		objective <- eval(substitute(mxFIMLObjective(covariance = x, means = y, thresholds = z),
			list(x = covName, y = meansName, z = .Object@thresholds)))
		model@objective <- objective
		return(model)
	}
)

mxRAMObjective <- function(A, S, F, M = NA, thresholds = NA) {
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
	if (is.na(thresholds)) thresholds <- as.integer(NA)
	return(new("MxRAMObjective", A, S, F, M, thresholds))
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
	if (is.na(objective@thresholds)) {
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

