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
		model <- listArgumentLISREL(model, lst, remove)
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
				stop("Unlike RAM, LISREL models must have latent variables.  The latentVars argument of a LISREL model must be a named list.")
			}
			if(length(names(manifestVars)) == 0 || !is.list(manifestVars) ){
				stop("The manifestVars argument of a LISREL model must be a named list.")
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
		theExpectation <- mxExpectationLISREL(LX='LX', PH='PH', TD='TD', TX=NA, KA=NA)
	# if has only endogenous variables
	} else if(length(mexog) == 0 & length(mendo) > 0){
		model <- addEndogenousMatrices(model, mendo, lendo)
		theExpectation <- mxExpectationLISREL(LY='LY', BE='BE', PS='PS', TE='TE', TY=NA, AL=NA)
	# if has both exogenous and endogenous variables
	} else {
		model <- addExoEndoMatrices(model, mexog, mendo, lexog, lendo)
		theExpectation <- mxExpectationLISREL(LX='LX', PH='PH', TD='TD', TX=NA, KA=NA, LY='LY', BE='BE', PS='PS', TE='TE', TY=NA, AL=NA, GA='GA', TH='TH')
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
				model <- addExogenousMeansLISREL(model)
			}
			if(hasEndogenousVars){
				model <- addEndogenousMeansLISREL(model)
			}
		}
	} else {
		model <- imxReplaceMethod(model, index, value)
	}
	return(model)
}

addExogenousMeansLISREL <- function(model){
	if (is.null(model[['TX']])) {
		model[['TX']] <- createMatrixLISREL(model, model@manifestVars$exogenous, "one", 'TX')
	}
	if (is.null(model[['KA']])) {
		model[['KA']] <- createMatrixLISREL(model, model@latentVars$exogenous, "one", 'KA')
	}
	model@expectation@TX <- "TX"
	model@expectation@KA <- "KA"
	return(model)
}

addEndogenousMeansLISREL <- function(model){
	if (is.null(model[['TY']])) {
		model[['TY']] <- createMatrixLISREL(model, model@manifestVars$endogenous, "one", 'TY')
	}
	if (is.null(model[['AL']])) {
		model[['AL']] <- createMatrixLISREL(model, model@latentVars$endogenous, "one", 'AL')
	}
	model@expectation@TY <- "TY"
	model@expectation@AL <- "AL"
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
	if (length(entries) == 0) {
		return(model)
	}
	
	filter <- sapply(entries, is, "MxPath")
	paths <- entries[filter]
	if (length(paths) > 0) {
		model <- insertAllPathsLISREL(model, paths)
	}
	filter <- sapply(entries, is, "MxThreshold")
	thresholds <- entries[filter]
	if(length(thresholds) > 0) {
		model <- insertAllThresholdsRAM(model, thresholds) # sic.  Re-use RAM threholds
	}
	filter <- sapply(entries, is, "MxData")
	data <- entries[filter]
	if (length(data) > 0) {
		if (length(data) > 1) {
			warning("Multiple data sources specified.  Only one will be chosen.")
		}
		data <- data[[1]]
		model@data <- data
		# If the data are WLS, then change the fit function to WLS away from the default ML.
		if(model@data@type=="acov" && class(model@fitfunction) %in% "MxFitFunctionML"){
			model[['fitfunction']] <- mxFitFunctionWLS()
		} else if(model@data@type %in% c('raw', 'cov') && !(class(model@fitfunction) %in% "MxFitFunctionML")){
			model[['fitfunction']] <- mxFitFunctionML()
		}
		#model[['F']] <- createMatrixF(model) #TODO something here about re-structuring LX and LY if needed
	}
	return(model)
}

insertAllPathsLISREL <-  function(model, paths){
	manivars <- model@manifestVars
	latevars <- model@latentVars
	exvars <- c(manivars$exogenous, latevars$exogenous)
	envars <- c(manivars$endogenous, latevars$endogenous)
	theVariables <- list(lex=latevars$exogenous, len=latevars$endogenous, mex=manivars$exogenous, men=manivars$endogenous)
	theMatrices <- list(LX=model[['LX']], LY=model[['LY']], BE=model[['BE']], GA=model[['GA']], PH=model[['PH']],
		PS=model[['PS']], TD=model[['TD']], TE=model[['TE']], TH=model[['TH']])
	TX <- model[['TX']]
	TY <- model[['TY']]
	KA <- model[['KA']]
	AL <- model[['AL']]
	
	legalVars <- c(exvars, envars, "one")
	
	for(i in 1:length(paths)) {
		path <- paths[[i]]
	
		missingvalues <- is.na(path@values)
		path@values[missingvalues] <- 0
		
		if (single.na(path@to)) {
			path@to <- path@from
			paths[[i]] <- path
		}
		
		allFromTo <- unique(c(path@from, path@to))
		varExist <- allFromTo %in% legalVars 
		if(!all(varExist)) {
			missingVars <- allFromTo[!varExist]
			stop(paste("Nice try, you need to add", 
				omxQuotes(missingVars), 
				"to either manifestVars or latentVars before you",
				"can use them in a path."), call. = FALSE)
		}
		
		if (length(path@from) == 1 && (path@from == "one")) {
			if ( (is.null(TX) || is.null(KA)) && any(path@to %in% exvars) ) {
				model <- addExogenousMeansLISREL(model)
				TX <- model[['TX']]
				KA <- model[['KA']]
			}
			if ( (is.null(TY) || is.null(AL)) && any(path@to %in% envars) ) {
				model <- addEndogenousMeansLISREL(model)
				TY <- model[['TY']]
				AL <- model[['AL']]
			}
			LISRELMeans <- insertMeansPathLISREL(path, TX, TY, KA, AL, manifest=manivars, latent=latevars)
			TX <- LISRELMeans[[1]]; TY <- LISRELMeans[[2]]; KA <- LISRELMeans[[3]]; AL <- LISRELMeans[[4]]
			model <- updateLISRELMeans(model, LISRELMeans)
		} else {
			expanded <- expandPathConnect(path@from, path@to, path@connect)
			path@from <- expanded$from
			path@to   <- expanded$to
			theMatrices <- insertPathLISREL(path, theMatrices, theVariables)
			LISRELMeans <- NULL
		}
	}
	
	model <- updateLISRELMatrices(model, theMatrices)
	
	return(model)
}

insertPathLISREL <- function(path, matrices, variables){
	allfrom <- path@from
	allto <- path@to
	allarrows <- path@arrows
	allfree <- path@free
	allvalues <- path@values
	alllabels <- path@labels
	alllbound <- path@lbound
	allubound <- path@ubound
	maxlength <- max(length(allfrom), length(allto))
	for(i in 0:(maxlength - 1)) {
		from <- allfrom[[i %% length(allfrom) + 1]]
		to <- allto[[i %% length(allto) + 1]]
		arrows <- allarrows[[i %% length(allarrows) + 1]]
		new <- list()
		new$value <- allvalues[[i %% length(allvalues) + 1]]
		new$free <- allfree[[i %% length(allfree) + 1]]
		new$label <- alllabels[[i %% length(alllabels) + 1]]
		new$ubound <- allubound[[i %% length(allubound) + 1]]
		new$lbound <- alllbound[[i %% length(alllbound) + 1]]
		
		#N.B. assuming that length(from) and length(to) are both 1
		if(from %in% variables$lex){
			matrices <- insertLatentExoPathsLISREL(from, to, arrows, old=matrices, new, variables)
		} else if(from %in% variables$len){
			matrices <- insertLatentEndoPathsLISREL(from, to, arrows, old=matrices, new, variables)
		} else if(from %in% variables$mex){
			matrices <- insertManifestExoPathsLISREL(from, to, arrows, old=matrices, new, variables)
		} else if(from %in% variables$men){
			matrices <- insertManifestEndoPathsLISREL(from, to, arrows, old=matrices, new, variables)
		}
	}
	return(matrices)
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

assignNewMatrixStuff <- function(from, to, oldMat, newStuff, remove=FALSE){
	if(remove==FALSE){
		oldMat$values[to, from] <- newStuff$value
		oldMat$free[to, from] <- newStuff$free
		oldMat$labels[to, from] <- newStuff$label
		oldMat$lbound[to, from] <- newStuff$lbound
		oldMat$ubound[to, from] <- newStuff$ubound
	} else {
		oldMat$values[to, from] <- 0
		oldMat$labels[to, from] <- as.character(NA)
		oldMat$free[to, from] <- FALSE
	}
	return(oldMat)
}

insertLatentExoPathsLISREL <- function(from, to, arrows, old, new, variables){
	if(to %in% variables$lex){
		if(arrows==1){
			stop('In LISREL, one-headed arrows do not exist between latent exogenous variables.  They are allowed when the destination is endogenous, and two-headed arrows are okay.')
		} else if(arrows==2){
			#add to PH
			old$PH <- assignNewMatrixStuff(from, to, old$PH, new)
			old$PH <- assignNewMatrixStuff(to, from, old$PH, new)
		}
	} else if(to %in% variables$len){
		if(arrows==1){
			#add to GA
			old$GA <- assignNewMatrixStuff(from, to, old$GA, new)
		} else if(arrows==2){
			stop('In LISREL, two-headed arrows do not exist between latent exogenous and latent endogenous variables.  One-headed arrows are okay.')
		}
	} else if(to %in% variables$mex){
		if(arrows==1){
			#add to LX
			old$LX <- assignNewMatrixStuff(from, to, old$LX, new)
		} else if(arrows==2){
			stop('In LISREL, two-headed arrows do not exist between latent and manifest exogenous variables.')
		}
	} else if(to %in% variables$men){
		stop('In LISREL, one-headed arrows cannot go from latent exogenous variables to manifest endogenous variables.')
	}
	return(old)
}

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

insertLatentEndoPathsLISREL <- function(from, to, arrows, old, new, variables){
	if(to %in% variables$lex){
		stop('In LISREL, no paths can go from latent endogenous variables to latent exogenous variables.')
	} else if(to %in% variables$len){
		if(arrows==1){
			#add to BE
			old$BE <- assignNewMatrixStuff(from, to, old$BE, new)
			old$PS <- assignNewMatrixStuff(from, to, old$PS, new, remove=TRUE)
		} else if(arrows==2){
			#add to PS
			old$PS <- assignNewMatrixStuff(from, to, old$PS, new)
			old$PS <- assignNewMatrixStuff(to, from, old$PS, new)
			old$BE <- assignNewMatrixStuff(from, to, old$BE, new, remove=TRUE)
		}
	} else if(to %in% variables$mex){
		stop('In LISREL, no paths can go from latent endogenous variables to manifest exogenous variables.')
	} else if(to %in% variables$men){
		if(arrows==1){
			#add to LY
			old$LY <- assignNewMatrixStuff(from, to, old$LY, new)
		} else if(arrows==2){
			stop('In LISREL, two-headed arrows do not exist between latent and manifest endogenous variables.')
		}
	}
	return(old)
}

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

insertManifestExoPathsLISREL <- function(from, to, arrows, old, new, variables){
	if(to %in% variables$lex){
		stop('In LISREL, paths are not allowed from manifest variables to latent variables.')
	} else if(to %in% variables$len){
		stop('In LISREL, paths are not allowed from manifest variables to latent variables.')
	} else if(to %in% variables$mex){
		if(arrows==1){
			stop('In LISREL, one-headed arrows are not allowed between manifest exogenous variables.')
		} else if(arrows==2){
			#add to TD
			old$TD <- assignNewMatrixStuff(from, to, old$TD, new)
			old$TD <- assignNewMatrixStuff(to, from, old$TD, new)
		}
	} else if(to %in% variables$men){
		if(arrows==1){
			stop('In LISREL, one-headed arrows are not allowed between manifest exogenous and manifest endogenous variables.')
		} else if(arrows==2){
			#add to TH
			old$TH <- assignNewMatrixStuff(from, to, old$TH, new)
			# N.B. TH is asymmetric so do NOT at the to->from path
		}
	}
	return(old)
}

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

insertManifestEndoPathsLISREL <- function(from, to, arrows, old, new, variables){
	if(to %in% variables$lex){
		stop('In LISREL, paths are not allowed from manifest variables to latent variables.')
	} else if(to %in% variables$len){
		stop('In LISREL, paths are not allowed from manifest variables to latent variables.')
	} else if(to %in% variables$mex){
		if(arrows==1){
			stop('In LISREL, one-headed arrows are not allowed between manifest exogenous variables.')
		} else if(arrows==2){
			#add to TH
			old$TH <- assignNewMatrixStuff(to, from, old$TH, new) #Note the switch of the from/to
			# N.B. TH is asymmetric so do NOT at the to->from path
		}
	} else if(to %in% variables$men){
		if(arrows==1){
			stop('In LISREL, one-headed arrows are not allowed between manifest variables.')
		} else if(arrows==2){
			#add to TE
			old$TE <- assignNewMatrixStuff(from, to, old$TE, new)
			old$TE <- assignNewMatrixStuff(to, from, old$TE, new)
		}
	}
	return(old)
}

# 'one' -> latentExo
#	if arrows==1
#		add to KA
#	if arrows==2
#		stop('nonsense')
# 'one' -> latentEndo
#	if arrows==1
#		add to AL
#	if arrows==2
#		stop('nonsense')
# 'one' -> manifestExo
#	if arrows==1
#		add to TX
#	if arrows==2
#		stop('nonsense')
# 'one' -> manifestEndo
#	if arrows==1
#		add to TY
#	if arrows==2
#		stop('nonsense')
# 'one' -> 'one'
#	stop('nonsense')

insertMeansPathLISREL <- function(path, TX, TY, KA, AL, manifest, latent){
	mexog <- manifest$exogenous
	mendo <- manifest$endogenous
	lexog <- latent$exogenous
	lendo <- latent$endogenous
	allto     <- path@to
	allarrows <- path@arrows
	allfree   <- path@free
	allvalues <- path@values
	alllabels <- path@labels
	alllbound <- path@lbound
	allubound <- path@ubound
	if (any(allarrows != 1)) {
		stop(paste('The means path must be a single-headed arrow\n',
		'path from "one" to variable(s)', omxQuotes(allto)), call. = FALSE)
	}
	for(i in 0:(length(allto) - 1)) {
		to <- allto[[i %% length(allto) + 1]]
		nextvalue  <- allvalues[[i %% length(allvalues) + 1]]
		nextfree   <- allfree[[i %% length(allfree) + 1]]
		nextlabel  <- alllabels[[i %% length(alllabels) + 1]]
		nextubound <- allubound[[i %% length(allubound) + 1]]
		nextlbound <- alllbound[[i %% length(alllbound) + 1]]
		if(to %in% mexog){ #Update TX
			TX@free[to, 1] <- nextfree
			TX@values[to, 1] <- nextvalue
			TX@labels[to, 1] <- nextlabel
			TX@ubound[to, 1] <- nextubound
			TX@lbound[to, 1] <- nextlbound
		} else if (to %in% mendo){ #Update TY
			TY@free[to, 1] <- nextfree
			TY@values[to, 1] <- nextvalue
			TY@labels[to, 1] <- nextlabel
			TY@ubound[to, 1] <- nextubound
			TY@lbound[to, 1] <- nextlbound
		} else if (to %in% lexog){ #Update KA
			KA@free[to, 1] <- nextfree
			KA@values[to, 1] <- nextvalue
			KA@labels[to, 1] <- nextlabel
			KA@ubound[to, 1] <- nextubound
			KA@lbound[to, 1] <- nextlbound
		} else if (to %in% lendo){ #Update AL
			AL@free[to, 1] <- nextfree
			AL@values[to, 1] <- nextvalue
			AL@labels[to, 1] <- nextlabel
			AL@ubound[to, 1] <- nextubound
			AL@lbound[to, 1] <- nextlbound
		} else {
			stop("I can't figure out where to put the means paths.")
		}
	}
	return(list(TX, TY, KA, AL))
}

updateLISRELMeans <- function(model, means){
	if(!single.na(means[[1]])){
		model[['TX']] <- means[[1]]
	}
	if(!single.na(means[[2]])){
		model[['TY']] <- means[[2]]
	}
	if(!single.na(means[[3]])){
		model[['KA']] <- means[[3]]
	}
	if(!single.na(means[[4]])){
		model[['AL']] <- means[[4]]
	}
	return(model)
}

updateLISRELMatrices <- function(model, matrices){
	if(!single.na(matrices$LX)){
		model[['LX']] <- matrices$LX
	}
	if(!single.na(matrices$LY)){
		model[['LY']] <- matrices$LY
	}
	if(!single.na(matrices$BE)){
		model[['BE']] <- matrices$BE
	}
	if(!single.na(matrices$GA)){
		model[['GA']] <- matrices$GA
	}
	if(!single.na(matrices$PH)){
		model[['PH']] <- matrices$PH
	}
	if(!single.na(matrices$PS)){
		model[['PS']] <- matrices$PS
	}
	if(!single.na(matrices$TD)){
		model[['TD']] <- matrices$TD
	}
	if(!single.na(matrices$TE)){
		model[['TE']] <- matrices$TE
	}
	if(!single.na(matrices$TH)){
		model[['TH']] <- matrices$TH
	}
	return(model)
}
