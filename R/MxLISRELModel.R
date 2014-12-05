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
		# Returns an ML fitfunction and an expectation with al NA matrices
		# add needed matrices later
		if (is.null(model[['expectation']])) {
			model[['expectation']] <- mxExpectationLISREL()
		}
		if (is.null(model[['fitfunction']])) {
			model[['fitfunction']] <- mxFitFunctionML()
		}
		return(model)
	}
)

addExogenousMatrices <- function(model, manifest, latent){
	if (is.null(model[['LX']])) {
		model[['LX']] <- createMatrixLISREL(model, manifest, latent, 'LX')
	}
	if (is.null(model[['PH']])) {
		model[['PH']] <- createMatrixLISREL(model, latent, latent, 'PH')
	}
	if (is.null(model[['TD']])) {
		model[['TD']] <- createMatrixLISREL(model, manifest, manifest, 'TD')
	}
	# Do we add means???
	if (is.null(model[['TX']])) {
		model[['TX']] <- createMatrixLISREL(model, latent, NA, 'TX')
	}
	if (is.null(model[['KA']])) {
		model[['KA']] <- createMatrixLISREL(model, manifest, NA, 'KA')
	}
	return(model)
}

addEndogenousMatrices <- function(model, manifest, latent){
	if (is.null(model[['LY']])) {
		model[['LY']] <- createMatrixLISREL(model, manifest, latent, 'LY')
	}
	if (is.null(model[['BE']])) {
		model[['BE']] <- createMatrixLISREL(model, latent, latent, 'BE')
	}
	if (is.null(model[['PS']])) {
		model[['PS']] <- createMatrixLISREL(model, latent, latent, 'PS')
	}
	if (is.null(model[['TE']])) {
		model[['TE']] <- createMatrixLISREL(model, manifest, manifest, 'TE')
	}
	# Do we add means???
	if (is.null(model[['TY']])) {
		model[['TY']] <- createMatrixLISREL(model, latent, NA, 'TY')
	}
	if (is.null(model[['AL']])) {
		model[['AL']] <- createMatrixLISREL(model, manifest, NA, 'AL')
	}
	return(model)
}

addExoEndoMatrices <- function(model, manifestExo, manifestEndo, latentExo, latentEndo){
	model <- addExogenousMatrices(model, manifestExo, latentExo)
	model <- addEndogenousMatrices(model, manifestEndo, latentEndo)
	if (is.null(model[['GA']])) {
		model[['GA']] <- createMatrixLISREL(model, latentEndo, latentExo, 'GA')
	}
	if (is.null(model[['TH']])) {
		model[['TH']] <- createMatrixLISREL(model, manifestEndo, manifestExo, 'TH')
	}
	return(model)
}

createMatrixLISREL <- function(model, rowvariables, colvariables, matrixname){
	clen <- length(colvariables)
	rlen <- length(rowvariables)
	names <- list(rowvariables, colvariables)
	values <- matrix(0, rlen, clen)
	free <- matrix(FALSE, rlen, clen)
	labels <- matrix(as.character(NA), rlen, clen)
	if(matrixname %in% c('LX', 'LY', 'BE', 'GA', 'TX', 'TY', 'KA', 'AL')){
		matrixtype <- 'Full'
	} else if(matrixname %in% c('PH', 'PS', 'TD', 'TE', 'TH')){	
		matrixtype <- 'Symm'
	}
	retval <- mxMatrix(matrixtype, values = values, free = free, labels = labels, name = matrixname)
	dimnames(retval) <- names
	return(retval)
}



setMethod("imxModelBuilder", "MxLISRELModel", 
	function(model, lst, name, manifestVars, latentVars, submodels, remove, independent) {
		model <- nameArgument(model, name)
		model <- variablesArgumentLISREL(model, manifestVars, latentVars, submodels, remove)
#		model <- listArgumentLISREL(model, lst, remove)
		notPathOrData <- getNotPathsOrData(lst)
		callNextMethod(model, notPathOrData, NA, character(), character(), list(), remove, independent)
		#stop("Not implemented")
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
	# works for initialization
	# TODO add processing for adding/removing variables and variables switching from exo to endo.
	model <- addVariablesHelper(model, "latentVars", latent)
	model <- addVariablesHelper(model, "manifestVars", manifest)
	mexog <- model@manifestVars$exogenous
	mendo <- model@manifestVars$endogenous
	lexog <- model@latentVars$exogenous
	lendo <- model@latentVars$endogenous
	# if has only exogenous variables
	if(length(mexog) > 0 & length(mendo) == 0){
		model <- addExogenousMatrices(model, mexog, lexog)
		theExpectation <- mxExpectationLISREL(LX='LX', PH='PH', TD='TD', TX='TX', KA='KA') #do we add means?
	# if has only endogenous variables
	} else if(length(mexog) == 0 & length(mendo) > 0){
		model <- addEndogenousMatrices(model, mendo, lendo)
		theExpectation <- mxExpectationLISREL(LY='LY', BE='BE', PS='PS', TE='TE', TY='TY', AL='AL') #do we add means?
	# if has both exogenous and endogenous variables
	} else {
		model <- addExoEndoMatrices(model, mexog, mendo, lexog, lendo)
		theExpectation <- mxExpectationLISREL(LX='LX', PH='PH', TD='TD', TX='TX', KA='KA', LY='LY', BE='BE', PS='PS', TE='TE', TY='TY', AL='AL', GA='GA', TH='TH')
	}
	model[['expectation']] <- theExpectation
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

listArgumentLISREL <- function(model, lst, remove){
	if(remove == TRUE) {
		model <- removeEntriesLISREL(model, lst)
	} else {
		model <- addEntriesLISREL(model, lst)
	}
	return(model)
}

removeEntriesLISREL <- function(model, entries){
	stop("Removing paths not yet implemented for LISREL")
}

addEntriesLISREL <- function(model, entries){
	stop("Adding paths not yet implemented for LISREL")
}

# Schematic of where each entry goes depending on its 'from' and 'to' and 'arrows' arguments
#from -> to
# We could first try to assign the path to one of the 13 valid combinations
# and then throw an error with the particulars of the invalid combination.
# stop(paste('Tried to add a', oneOrtwo, '-headed path from an', exoORendoFrom, latentORmanifestFrom, 'variable to an',
#	exoORendoTo, latentORmanifestTo, 'variable, but that is not possible in LISREL.'))
#latentExo -> latentExo
#	if arrows==1
#		stop('nonsense')
#	if arrows==2
#		add to PH
#latentExo -> latentEndo
#	if arrows==1
#		add to GA
#	if arrows=2
#		stop('nonsense')
#latentExo -> manifestExo
#	if arrows=1
#		add to LX
#	if arrows=2
#		stop('nonsense')
#latentExo -> manifestEndo
#	stop('nonsense')
#latentExo -> 'one'
#	stop('nonsense')

#latentEndo -> latentExo
#	stop('nonsense')
#latentEndo -> latentEndo
#	if arrows==1
#		add to BE
#	if arrows==2
#		add to PS
#latentEndo -> manifestExo
#	stop('nonsense')
#latentEndo -> manifestEndo
#	if arrows==1
#		add to LY
#	if arrows==2
#		stop('nonsense')
#latentEndo -> 'one'
#	stop('nonsense')

#manifestExo -> latentExo
#	stop('nonsense')
#manifestExo -> latentEndo
#	stop('nonsense')
#manifestExo -> manifestExo
#	if arrows==1
#		stop('nonsense')
#	if arrows==2
#		add to TD
#manifestExo -> manifestEndo
#	if arrows==1
#		stop('nonsense')
#	if arrows==2
#		add to TH
#manifestExo -> 'one'
#	stop('nonsense')

#manifestEndo -> latentExo
#	stop('nonsense')
#manifestEndo -> latentEndo
#	stop('nonsense')
#manifestEndo -> manifestExo
#	if arrows==1
#		stop('nonsense')
#	if arrows==2
#		add to TH
#manifestEndo -> manifestEndo
#	if arrows==1
#		stop('nonsense')
#	if arrows==2
#		add to TE
#manifestEndo -> 'one'
#	stop('nonsense')

#'one' -> latentExo
#	if arrows==1
#		add to KA
#	if arrows==2
#		stop('nonsense')
#'one' -> latentEndo
#	if arrows==1
#		add to AL
#	if arrows==2
#		stop('nonsense')
#'one' -> manifestExo
#	if arrows==1
#		add to TX
#	if arrows==2
#		stop('nonsense')
#'one' -> manifestEndo
#	if arrows==1
#		add to TY
#	if arrows==2
#		stop('nonsense')
#'one' -> 'one'
#	stop('nonsense')


