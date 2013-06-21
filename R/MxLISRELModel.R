#
#   Copyright 2007-2013 The OpenMx Project
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

setClass(Class = "MxLISRELModel",
	representation = representation(),
	contains = "MxModel")

imxModelTypes[['LISREL']] <- "MxLISRELModel"

# imxVariableTypes <- c(imxVariableTypes, "exogenous", "endogenous")


#--------------------------------------------------------------------
# Define generic functions

setMethod("imxTypeName", "MxLISRELModel", 
	function(model) { "LISREL" }
)

setMethod("imxInitModel", "MxLISRELModel", 
	function(model) {
		stop("Not implemented")
		# Returns an ML fitfunction and an expectation with all NA matrices
		# Then later add matrices when I know what dims they have.
		if (is.null(model[['expectation']])) {
			model[['expectation']] <- mxExpectationLISREL()
		}
		if (is.null(model[['fitfunction']])) {
			model[['fitfunction']] <- mxFitFunctionML()
		}
		return(model)
	}
)


createMatrixLX <- function(model){
	lvariables <- c(model@latentVars$exogenous)
	mvariables <- c(model@maifestVars$exogenous)
	llen <- length(lvariables)
	mlen <- length(mvariables)
	names <- list(mvariables, lvariables)
	values <- matrix(0, mlen, llen)
	free <- matrix(FALSE, mlen, llen)
	labels <- matrix(as.character(NA), mlen, llen)
	retval <- mxMatrix("Full", values = values, free = free, labels = labels, name = "LX")
	dimnames(retval) <- names
	return(retval)
}

createMatrixLY <- function(model){
	lvariables <- c(model@latentVars$endogenous)
	mvariables <- c(model@maifestVars$endogenous)
	llen <- length(lvariables)
	mlen <- length(mvariables)
	names <- list(mvariables, lvariables)
	values <- matrix(0, mlen, llen)
	free <- matrix(FALSE, mlen, llen)
	labels <- matrix(as.character(NA), mlen, llen)
	retval <- mxMatrix("Full", values = values, free = free, labels = labels, name = "LY")
	dimnames(retval) <- names
	return(retval)
}

createMatrixBE <- function(model){
	variables <- c(model@latentVars$endogenous)
	len <- length(variables)
	names <- list(variables, variables)
	values <- matrix(0, len, len)
	free <- matrix(FALSE, len, len)
	labels <- matrix(as.character(NA), len, len)
	retval <- mxMatrix("Full", values = values, free = free, labels = labels, name = "BE")
	dimnames(retval) <- names
	return(retval)
}

createMatrixGA <- function(model){
	xvariables <- c(model@latentVars$exogenous)
	yvariables <- c(model@latentVars$endogenous)
	xlen <- length(xvariables)
	ylen <- length(yvariables)
	names <- list(yvariables, xvariables)
	values <- matrix(0, ylen, xlen)
	free <- matrix(FALSE, ylen, xlen)
	labels <- matrix(as.character(NA), ylen, xlen)
	retval <- mxMatrix("Full", values = values, free = free, labels = labels, name = "GA")
	dimnames(retval) <- names
	return(retval)
}

# TODO Fill in the rest of the createMatrix* functions

createMatrixPH <- function(model){} #Latent cov of xi
createMatrixPS <- function(model){} #Latent cov of eta
createMatrixTD <- function(model){} #residu cov of x
createMatrixTE <- function(model){} #residu cov of y
createMatrixTH <- function(model){} #residu cov of xy
createMatrixTX <- function(model){} #interc of x
createMatrixTY <- function(model){} #interc of y
createMatrixKA <- function(model){} #mean of xi
createMatrixAL <- function(model){} #mean of eta




setMethod("imxModelBuilder", "MxLISRELModel", 
	function(model, lst, name, 
		 manifestVars, latentVars, submodels, remove, independent) {
		model <- nameArgument(model, name)
		model <- variablesArgumentLISREL()
		model <- listArgumentLISREL(model, lst, remove)
		notPathOrData <- getNotPathsOrData(lst)
		callNextMethod()
		stop("Not implemented")
	}
)

setMethod("imxVerifyModel", "MxLISRELModel",
	function(model) {
		# TODO somewhere in here add check that at least one of LX or LY exist
		return(TRUE)
	}
)


setReplaceMethod("[[", "MxLISRELModel",
	function(x, i, j, value) {
		stop("Not implemented")
	}
)

setReplaceMethod("$", "MxLISRELModel",
	function(x, name, value) {
		stop("Not implemented")
	}
)

# Helpers for LISREL models
variablesArgumentLISREL <- function(model, manifestVars, latentVars, submodels, remove){
	# include check prior to varsToCharacter to see if length(names(latentVars)) == 0 (and if same for manifestVars
}
