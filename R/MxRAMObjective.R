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
		M = "MxCharOrNumber"),
	contains = "MxBaseObjective")

setMethod("initialize", "MxRAMObjective",
	function(.Object, A, S, F, M,  
		data = as.numeric(NA), name = 'objective') {
		.Object@name <- name
		.Object@A <- A
		.Object@S <- S
		.Object@F <- F
		.Object@M <- M
		.Object@data <- data
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
			msg <- paste("The MxRAMObjective", omxQuotes(name),
				"does not have a dataset associated with it in model",
				omxQuotes(flatModel@name))
			stop(msg, call. = FALSE)
		}
		.Object@A <- omxLocateIndex(flatModel, aMatrix, name)
		.Object@S <- omxLocateIndex(flatModel, sMatrix, name)
		.Object@F <- omxLocateIndex(flatModel, fMatrix, name)
		if (is.na(.Object@M)) {
			.Object@M <- as.integer(NA)
		} else {
			.Object@M <- omxLocateIndex(flatModel, mMatrix, name)
		}
		.Object@data <- omxLocateIndex(flatModel, data, name)
		return(.Object)
})

setMethod("omxObjModelConvert", "MxRAMObjective",
	function(.Object, model) {
		if (is.null(model[['data']])) {
			return(NA)
		} else if (model[['data']]@type != 'raw') {
			return(NA)
		} else if (is.na(.Object@M)) {
			return(NA)
		}
		if (is.null(model[['I']])) {
			iName <- 'I'
		} else {
			iName <- omxUntitledName()
		}
		model <- mxModel(model, mxMatrix(type="Iden", 
			nrow=nrow(model[[.Object@A]]), name = iName))
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
		manifestVars <- dimnames(model[['F']])[[1]]
		dimnames(algebra) <- list(manifestVars, manifestVars)
		model <- mxModel(model, algebra)	
		objective <- eval(substitute(mxFIMLObjective(x, y),
			list(x = covName, y = .Object@M)))
		model@objective <- objective
		return(model)
	}
)

mxRAMObjective <- function(aMatrix = "A", sMatrix = "S", fMatrix = "F", mMatrix = NA) {
	if (typeof(aMatrix) != "character") {
		msg <- paste("aMatrix argument is not a string",
			"(the name of the 'A' matrix)")
		stop(msg)
	}	
	if (typeof(sMatrix) != "character") {
		msg <- paste("sMatrix argument is not a string",
			"(the name of the 'S' matrix)")
		stop(msg)
	}
	if (typeof(fMatrix) != "character") {
		msg <- paste("fMatrix argument is not a string",
			"(the name of the 'F' matrix)")
		stop(msg)
	}
	if (is.na(mMatrix)) mMatrix <- as.character(NA)
	if (typeof(mMatrix) != "character") {
		msg <- paste("mMatrix argument is not a string",
			"(the name of the 'M' matrix)")
		stop(msg)
	}	
	return(new("MxRAMObjective", aMatrix, sMatrix, fMatrix, mMatrix))
}

displayRAMObjective <- function(objective) {
	cat("MxRAMObjective", omxQuotes(objective@name), '\n')
	cat("A matrix :", omxQuotes(objective@A), '\n')
	cat("S matrix :", omxQuotes(objective@S), '\n')
	cat("F matrix :", omxQuotes(objective@F), '\n')
	if (is.na(objective@M)) {
		cat("M matrix :", objective@M, '\n')
	} else {
		cat("M matrix :", omxQuotes(objective@M), '\n')
	}
	if (length(objective@result) == 0) {
		cat("Result : empty\n")
	} else {
		cat("Result : \n") 
		print(objective@result)
	}
	invisible(objective)
}

setMethod("print", "MxRAMObjective", function(x,...) { 
	displayRAMObjective(x) 
})

setMethod("show", "MxRAMObjective", function(object) { 
	displayRAMObjective(object) 
})

