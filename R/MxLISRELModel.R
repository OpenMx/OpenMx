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
		
		if (is.null(model[['LX']])) {
			model[['LX']] <- createMatrixLX(model)
		}
		if (is.null(model[['LY']])) {
			model[['LY']] <- createMatrixLY(model)
		}
		if (is.null(model[['BE']])) {
			model[['BE']] <- createMatrixBE(model)
		}
		if (is.null(model[['GA']])) {
			model[['GA']] <- createMatrixGA(model)
		}
		if (is.null(model[['PH']])) {
			model[['PH']] <- createMatrixPH(model)
		}
		if (is.null(model[['PS']])) {
			model[['PS']] <- createMatrixPS(model)
		}
		if (is.null(model[['TD']])) {
			model[['TD']] <- createMatrixTD(model)
		}
		if (is.null(model[['TE']])) {
			model[['TE']] <- createMatrixTE(model)
		}
		if (is.null(model[['TH']])) {
			model[['TH']] <- createMatrixTH(model)
		}
		if (is.null(model[['TX']])) {
			model[['TX']] <- createMatrixTX(model)
		}
		if (is.null(model[['TY']])) {
			model[['TY']] <- createMatrixTY(model)
		}
		if (is.null(model[['KA']])) {
			model[['KA']] <- createMatrixKA(model)
		}
		if (is.null(model[['AL']])) {
			model[['AL']] <- createMatrixAL(model)
		}
		if (is.null(model[['expectation']])) {
			model[['expectation']] <- mxExpectationLISREL('LX', 'LY', 'BE', 'GA', 'PH', 'PS', 'TD', 'TE', 'TH', 'TX', 'TY', 'KA', 'AL')
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

# TODO See if there is a way to change an mxMatrix's type.  E.g. TD & TE are often but not always diagonal, and should(?) be stored as diagonal if possible.


setMethod("imxModelBuilder", "MxLISRELModel", 
	function(model, lst, name, 
		manifestVars, latentVars, submodels, remove, independent) {
		stop("Not implemented")
	}
)

setMethod("imxVerifyModel", "MxLISRELModel",
	function(model) {
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
