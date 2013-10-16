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

##' MxLISRELModel
##'
##' This is an internal class and should not be used directly.
##'
##' @aliases
##' $<-,MxLISRELModel-method
##' [[<-,MxLISRELModel-method
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
	function(model, lst, name, manifestVars, latentVars, submodels, remove, independent) {
		model <- nameArgument(model, name)
		model <- variablesArgumentLISREL(model, manifestVars, latentVars, submodels, remove)
		model <- listArgumentLISREL(model, lst, remove)
		notPathOrData <- getNotPathsOrData(lst)
		callNextMethod(model, notPathOrData, NA, character(), character(), list(), remove, independent)
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
		return(replaceMethodLISREL(x, i, value))
	}
)

setReplaceMethod("$", "MxLISRELModel",
	function(x, name, value) {
		stop("Not implemented")
	}
)

# Helpers for LISREL models
variablesArgumentLISREL <- function(model, manifestVars, latentVars, submodels, remove){
	if (single.na(manifestVars)) {
		manifestVars <- character()
	}
	if (single.na(latentVars)) {
		latentVars <- character()
	}
	if (remove == TRUE) {
		if (length(latentVars) + length(manifestVars) > 0) {
			model <- removeVariablesLISREL(model, latentVars, manifestVars)
		}
		if (length(submodels)) for(i in 1:length(submodels)) {
			model <- removeSingleNamedEntity(model, submodels[[i]])
		}
	} else {
		if (length(manifestVars) + length(latentVars) > 0) {
			if(length(names(latentVars)) == 0 || !is.list(latentVars) ){
				stop("The latentVars argument of a LISREL models must be a named list.")
			}
			if(length(names(manifestVars)) == 0 || !is.list(manifestVars) ){
				stop("The manifestVars argument of a LISREL models must be a named list.")
			}
			latentVars <- varsToCharacter(latentVars, "latent")
			manifestVars <- varsToCharacter(manifestVars, "manifest")
			checkVariables(model, latentVars, manifestVars)
			model <- addVariablesLISREL(model, latentVars, manifestVars)
		}
		if (length(submodels)) for(i in 1:length(submodels)) {
			model <- addSingleNamedEntity(model, submodels[[i]])
		}
	}
	return(model)
	# include check prior to varsToCharacter to see if length(names(latentVars)) == 0 (and if same for manifestVars
}

removeVariablesLISREL <- function(model, latent, manifest){
	return(model)
}

addVariablesLISREL <- function(model, latent, manifest){
	return(model)
}

addVariablesMatrixLISREL <- function(oldMatrix, value, model, numNewRows, numNewCols, varType, varSubtype){
	slot(model, varType)[[varSubtype]] #e.g. get the manifestVars$exogenous with varType="manifestVars" and varSubtype="exogenous"
}

replaceMethodLISREL <- function(model, index, value){
	pair <- imxReverseIdentifier(model, index)
	namespace <- pair[[1]]
	name <- pair[[2]]
	if (namespace == model@name && name == "data") {
		hasExogenousVars <- length(c(model@manifestVars$exogenous, model@latentVars$exogenous)) > 0
		hasEndogenousVars <- length(c(model@manifestVars$endogenous, model@latentVars$endogenous)) > 0
		model@data <- value
		if (requireMeansVector(value)) {
			if(hasExogenousVars){
				model@expectation@TX <- "TX"
				model@expectation@KA <- "KA"
			}
			if(hasExogenousVars){
				model@expectation@TY <- "TY"
				model@expectation@AL <- "AL"
			}
		}
	} else {
		model <- imxReplaceMethod(model, index, value)
	}
	return(model)
}
