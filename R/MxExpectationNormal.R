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

setClass(Class = "BaseExpectationNormal",
	 contains = "MxBaseExpectation")

setClass(Class = "MxExpectationNormal",
	representation = representation(
		covariance = "MxCharOrNumber",
		means = "MxCharOrNumber",
		thresholds = "MxCharOrNumber",
		dims = "character",
		dataColumns = "numeric",
		thresholdColumns = "numeric",
		thresholdLevels = "numeric",
		threshnames = "character",
		ExpCov = "matrix",
		ExpMean = "matrix",
	        numStats = "numeric"),
	contains = "BaseExpectationNormal")

setMethod("initialize", "MxExpectationNormal",
	function(.Object, covariance, means, dims, thresholds, threshnames,
		data = as.integer(NA), name = 'expectation') {
		.Object@name <- name
		.Object@covariance <- covariance
		.Object@means <- means
		.Object@data <- data
		.Object@thresholds <- thresholds
		.Object@dims <- dims
		.Object@threshnames <- threshnames
		.Object@ExpCov <- matrix()
		.Object@ExpMean <- matrix()
		return(.Object)
	}
)

setMethod("qualifyNames", signature("MxExpectationNormal"), 
	  function(.Object, modelname, namespace) {
		.Object@name <- imxIdentifier(modelname, .Object@name)
		.Object@covariance <- imxConvertIdentifier(.Object@covariance, 
			modelname, namespace)
		.Object@means <- imxConvertIdentifier(.Object@means, 
			modelname, namespace)
		.Object@data <- imxConvertIdentifier(.Object@data, 
			modelname, namespace)
		.Object@thresholds <- sapply(.Object@thresholds,
			imxConvertIdentifier, modelname, namespace)
		return(.Object)
})

setMethod("genericExpDependencies", signature("MxExpectationNormal"),
	function(.Object, dependencies) {
	sources <- c(.Object@covariance, .Object@means, .Object@thresholds)
	sources <- sources[!is.na(sources)]
	dependencies <- imxAddDependency(sources, .Object@name, dependencies)
	return(dependencies)
})

setMethod("genericExpConvertEntities", "MxExpectationNormal",
	function(.Object, flatModel, namespace, labelsData) {
		flatModel <- updateExpectationDimnames(.Object, flatModel, labelsData)
		flatModel <- updateThresholdDimnames(.Object, flatModel, labelsData)
		return(flatModel)
	}
)

setMethod("genericExpRename", signature("MxExpectationNormal"),
	function(.Object, oldname, newname) {
		.Object@means <- renameReference(.Object@means, oldname, newname)
		.Object@covariance <- renameReference(.Object@covariance, oldname, newname)
		.Object@data <- renameReference(.Object@data, oldname, newname)
		.Object@thresholds <- sapply(.Object@thresholds, renameReference, oldname, newname)		
		return(.Object)
})

setMethod("genericExpGetPrecision", "MxExpectationNormal",
	function(.Object) {
		if(!single.na(.Object@thresholds)) {
			return(list(stepSize=0.1, iterations=3L, functionPrecision=1e-9))
		} else {
			callNextMethod();
		}
})

setMethod("genericGetExpected", signature("MxExpectationNormal"),
	function(.Object, model, what, defvar.row=1) {
		ret <- list()
		if ('covariance' %in% what) {
			covname <- .Object@covariance
			cov <- mxEvalByName(covname, model, compute=TRUE, defvar.row=defvar.row)
			dnames <- .Object$dims
			if(!single.na(dnames)){
				colnames(cov) <- dnames
				rownames(cov) <- dnames
			}
			ret[['covariance']] <- cov
		}
		if ('means' %in% what) {
			meanname <- .Object@means
			if(!single.na(meanname)){
				mean <- mxEvalByName(meanname, model, compute=TRUE, defvar.row=defvar.row)
				dnames <- .Object$dims
				if(!single.na(dnames)){
					colnames(mean) <- dnames
				}
			} else {mean <- matrix( , 0, 0)}
			ret[['means']] <- mean
		}
		if ('thresholds' %in% what) {
			thrname <- .Object@thresholds
			if(!single.na(thrname)){
				thr <- mxEvalByName(thrname, model, compute=TRUE, defvar.row=defvar.row)
				tnames <- .Object$threshnames
				if(!single.na(tnames)){
					colnames(thr) <- tnames
				}
			} else {thr <- matrix( , 0 , 0)}
			ret[['thresholds']] <- thr
		}
		ret
	})

setMethod("genericGetExpectedVector", signature("BaseExpectationNormal"),
	function(.Object, model, defvar.row=1) {
		ret <- genericGetExpected(.Object, model, c('covariance', 'means', 'thresholds'), defvar.row)
		cov <- ret[['covariance']]
		mns <- ret[['means']]
		if (is.null(mns)) stop("mns is null")
		thr <- ret[['thresholds']]
		if (is.null(thr)) stop("thresholds is null")
		v <- c(vech(cov), mns[!is.na(mns)], thr[!is.na(thr)])
		return(v)
})

imxGetExpectationComponent <- function(model, component, defvar.row=1)
{
	if(is.null(model$expectation) && (class(model$fitfunction) %in% "MxFitFunctionMultigroup") ){
		submNames <- sapply(strsplit(model$fitfunction$groups, ".", fixed=TRUE), "[", 1)
		got <- list()
		for(amod in submNames){
			got[[amod]] <- imxGetExpectationComponent(model[[amod]], component, defvar.row=1)
		}
		if(component=='vector'){got <- unlist(got)}
		got
	} else if (length(component) == 1 && component == 'vector') {
		genericGetExpectedVector(model$expectation, model, defvar.row)
	} else {
		got <- genericGetExpected(model$expectation, model, component, defvar.row)
		if (length(got) == 1) {
			got[[1]]
		} else {
			got
		}
	}
}

mxGetExpected <- imxGetExpectationComponent

sse <- function(x){sum(x^2)}

#' Estimate the Jacobian of manifest model with respect to parameters
#'
#' The manifest model excludes any latent variables or processes. For
#' RAM and LISREL models, the manifest model contains only the
#' manifest variables with free means, covariance, and thresholds.
#'
#' The Jacobian is estimated by the central finite difference.
#'
#' @param model an mxModel
#' @param defvar.row which row to use for definition variables
#' @return a matrix with manifests in the rows and original parameters in the columns
#' @seealso \link{mxGetExpected}
omxManifestModelByParameterJacobian <- function(model, defvar.row=1) {
	theParams <- omxGetParameters(model)
	numDeriv::jacobian(func=.mat2param, x=theParams, method.args=list(r=2), model=model, defvar.row=defvar.row)
}

mxCheckIdentification <- function(model, details=TRUE){
	notAllowedFits <- c("MxFitFunctionAlgebra", "MxFitFunctionRow", "MxFitFunctionR")
	if( class(model$fitfunction) %in% notAllowedFits ){
		msg <- paste("Identification check is not possible for models with", omxQuotes(notAllowedFits), 'fit functions.\n', "If you have a multigroup model, use mxFitFunctionMultigroup.")
		stop(msg, call.=FALSE)
	}
	eps <- 1e-17
	theParams <- omxGetParameters(model)
	jac <- omxManifestModelByParameterJacobian(model)
	# Check that rank of jac == length(theParams)
	rank <- qr(jac)$rank
	if(rank == length(theParams)){
		message("Model is locally identified")
		stat <- TRUE
	} else {
		message("Model is not locally identified")
		stat <- FALSE
	}
	if(details == TRUE){
		jacOC <- Null(t(jac)) # Orthogonal complement of t(jac), i.e. the basis for the null space of the column space of jac
		nidp <- names(theParams)[apply(jacOC, 1, sse) > eps] # non-identified free params have non-zero rows
		if(length(nidp) == 0) {
			nidp <- "None"
		}
	} else {
		nidp <- "Not Requested"
	}
	return(list(status=stat, jacobian=jac, non_identified_parameters=nidp))
}

.mat2param <- function(x, model, defvar.row=1){
  paramNames <- names(omxGetParameters(model))
  model <- omxSetParameters(model, values=x, labels=paramNames, free=TRUE)
  got <- mxGetExpected(model, 'vector', defvar.row)
  got
}

setGeneric("genericGenerateData",
	function(.Object, model, nrows) {
	return(standardGeneric("genericGenerateData"))
})

setMethod("genericGenerateData", signature("MxExpectationNormal"),
	function(.Object, model, nrows) {
		return(generateNormalData(model, nrows))
})

generateNormalData <- function(model, nrows){
	# Check for definition variables
	if(imxHasDefinitionVariable(model)){
		if(nrows == nrow(model@data@observed)){
			# Generate data row by row
			theCov <- imxGetExpectationComponent(model, "covariance")
			data <- matrix(NA, nrow=nrows, ncol=ncol(theCov))
			colnames(data) <- colnames(theCov)
			data <- as.data.frame(data)
			for(i in 1:nrows){
				theMeans <- imxGetExpectationComponent(model, "means", defvar.row=i)
				theCov <- imxGetExpectationComponent(model, "covariance", defvar.row=i)
				theThresh <- imxGetExpectationComponent(model, "thresholds", defvar.row=i)
				data[i,] <- mvtnorm::rmvnorm(1, theMeans, theCov)
				data[i,] <- ordinalizeDataHelper(data[i,], theThresh)
			}
			data <- ordinalizeDataHelper(data, theThresh, cut=FALSE)
		} else{
			stop("Definition variable(s) found, but the number of rows in the data do not match the number of rows requested for data generation.")
		}
	} else{
		#use generic functions and mvtnorm::rmvnorm() to generate data
		theMeans <- imxGetExpectationComponent(model, "means")
		theCov <- imxGetExpectationComponent(model, "covariance")
		theThresh <- imxGetExpectationComponent(model, "thresholds")
		data <- mvtnorm::rmvnorm(nrows, theMeans, theCov)
		colnames(data) <- colnames(theCov)
		data <- as.data.frame(data)
		data <- ordinalizeDataHelper(data, theThresh)
	}
	return(data)
}

ordinalizeDataHelper <- function(data, thresh, cut=TRUE){
	if( prod(dim(thresh)) != 0){
		ordvars <- colnames(thresh)
		for(avar in ordvars){
			delthr <- thresh[,avar]
			usethr <- !is.na(delthr)
			delthr <- delthr[usethr]
			if(!is.null(rownames(thresh))){
				levthr <- rownames(thresh)[usethr]
			} else {
				levthr <- 1:(sum(usethr)+1)
			}
			if(cut==TRUE){
				delvar <- cut(as.vector(data[,avar]), c(-Inf, delthr, Inf), labels=levthr)
			} else{
				delvar <- data[,avar]
			}
			data[,avar] <- mxFactor(delvar, levels=levthr)
		}
	}
	return(data)
}

generateRelationalData <- function(model, returnModel) {
	plan <- mxComputeSequence(list(
	    mxComputeOnce('expectation', 'distribution', 'flat'),
	    mxComputeReportExpectation()
	))

	modelE <- mxModel(model, plan)
	modelE$expectation$.rampart <- 0L
	modelE <- mxRun(modelE, silent=TRUE)
	dataEnv <- new.env()
	for (dName in names(modelE@runstate$datalist)) {
		modelName <- substr(dName, 1, nchar(dName)-5)  # remove .data
		assign(modelName, modelE@runstate$datalist[[ dName ]]@observed, envir=dataEnv)
	}
	ed <- modelE$expectation$debug
	layout <- ed$layout
	fmt <- paste0('g%0', ceiling(log10(ed$numGroups)), 'd')
	for (gx in 1:ed$numGroups) {
		groupName <- sprintf(fmt, gx)
		numCopies <- length(unique(layout[layout$group == gx, 'copy']))
		cxLength <- length(ed[[groupName]]$mean) / numCopies
		groupTodo <- ed$layout[ed[[groupName]]$layout[,'aIndex'],]
		for (cx in 1:numCopies) {
			todo <- groupTodo[groupTodo$copy == cx,]
			repl1 <- mvtnorm::rmvnorm(1, ed[[groupName]]$mean[seq(1+(cx-1)*cxLength, cx*cxLength)],
						  sigma=as.matrix(ed[[groupName]]$covariance))
			dx <- 1
			for (tx in 1:nrow(todo)) {
				modelName <- as.character(todo[tx,'model'])
				if (modelName == modelE$name) {
					submodel <- modelE
				} else {
					submodel <- modelE[[modelName]]
				}
				row <- todo[tx,'row']
				manifests <- rownames(submodel$F)
				beforeData <- dataEnv[[modelName]][row, manifests]
				notMissing <- !is.na(beforeData)
				if (sum(notMissing) > 0) {
					afterData <- repl1[seq(dx, dx+sum(notMissing) - 1)]
					dataEnv[[modelName]][row, manifests[notMissing] ] <- afterData
					dx <- dx + sum(notMissing)
				}
			}
		}
	}
	if (!returnModel) {
		ret <- list()
		for (n in names(dataEnv)) {
			ret[[n]] <- dataEnv[[n]]
		}
		ret
	} else {
		for (modelName in names(dataEnv)) {
			if (modelName == model$name) {
				model@data@observed <- dataEnv[[modelName]]
			} else {
				model[[modelName]]@data@observed <- dataEnv[[modelName]]
			}
		}
		model
	}
}

mxGenerateData <- function(model, nrows=NULL, returnModel=FALSE) {
	fellner <- is(model$expectation, "MxExpectationRAM") && length(model$expectation$between);
	if (!fellner) {
		if (missing(nrows)) nrows <- nrow(model@data@observed)
		data <- genericGenerateData(model$expectation, model, nrows)
		if (returnModel) {
			model@data@observed <- as.data.frame(data)
			model
		} else {
			as.data.frame(data)
		}
	} else {
		if (!missing(nrows)) {
			stop("Specification of the number of rows is not supported for relational models")
		}
		generateRelationalData(model, returnModel)
	}
}

verifyExpectedObservedNames <- function(data, covName, flatModel, modelname, objectiveName) {
	covariance <- flatModel[[covName]]
	if (is(covariance, "MxMatrix") && !identical(dim(covariance), dim(data))) {
		msg <- paste("The dimensions for the expected covariance matrix",
			"and the observed covariance matrix",
			"in the", objectiveName, "expectation function in model",
			omxQuotes(modelname), "are not identical.")
		stop(msg, call. = FALSE)
	}
	if (!identical(dimnames(covariance), dimnames(data))) {
		msg <- paste("The dimnames for the expected covariance matrix",
			"and the observed covariance matrix",
			"in the", objectiveName, "expectation function in model",
			omxQuotes(modelname), "are not identical.")
		stop(msg, call. = FALSE)		
	}
}

verifyMeans <- function(meansName, mxDataObject, flatModel, modelname) {
	means <- flatModel[[meansName]]
	if (!is.null(means)) {
		if(any(is.na(mxDataObject@means))) {
			msg <- paste("In model", omxQuotes(modelname),
				"the Normal expectation function contains an expected means",
				"vector but the model is missing some data",
				"for the observed means.")
			stop(msg, call. = FALSE)
		}
		meanDimnames <- dimnames(means)		
	}
}

setMethod("genericExpFunConvert", "MxExpectationNormal", 
	function(.Object, flatModel, model, labelsData, dependencies) {
		modelname <- imxReverseIdentifier(model, .Object@name)[[1]]
		name <- .Object@name
		if(is.na(.Object@data)) {
			msg <- paste("The normal expectation function",
				"does not have a dataset associated with it in model",
				omxQuotes(modelname))
			stop(msg, call.=FALSE)
		}

		mxDataObject <- flatModel@datasets[[.Object@data]]
		dataName <- .Object@data
		.Object@data <- imxLocateIndex(flatModel, .Object@data, name)
		threshName <- .Object@thresholds
		.Object@thresholds <- imxLocateIndex(flatModel, threshName, name)
		covName <- .Object@covariance
		covariance <- flatModel[[covName]]
		.Object@covariance <- imxLocateIndex(flatModel, .Object@covariance, name)
		meansName <- .Object@means
		.Object@means <- imxLocateIndex(flatModel, .Object@means, name)

		if (inherits(mxDataObject, "MxDataDynamic")) return(.Object)

		if (mxDataObject@type != "raw") {
			verifyExpectedObservedNames(mxDataObject@observed, covName, flatModel, modelname, "Normal")
			verifyMeans(meansName, mxDataObject, flatModel, modelname)
		}
		verifyObservedNames(mxDataObject@observed, mxDataObject@means, mxDataObject@type, flatModel, modelname, "Normal")
		checkNumericData(mxDataObject)
		checkNumberOrdinalColumns(mxDataObject)
		covNames <- colnames(covariance)
		verifyMvnNames(covName, meansName, "expected", flatModel, modelname, class(.Object))
		.Object@dataColumns <- generateDataColumns(flatModel, covNames, dataName)
		verifyThresholds(flatModel, model, labelsData, dataName, covNames, threshName)
		retval <- generateThresholdColumns(flatModel, model, labelsData, covNames, dataName, threshName)
		.Object@thresholdColumns <- retval[[1]] 
		.Object@thresholdLevels <- retval[[2]]
		if (single.na(.Object@dims)) {
			.Object@dims <- covNames
		}
		return(.Object)
})

verifyMvnNames <- function(covName, meansName, type, flatModel, modelname, expectationName) {
	if (is.na(meansName)) {
		means <- NA
	} else {
		means <- flatModel[[meansName]]
	}
	covariance <- flatModel[[covName]]
	if (length(covariance)) {
		covDimnames <- dimnames(covariance)
		if (is.null(covDimnames)) {
			msg <- paste("The",type,"covariance matrix associated",
				     "with", expectationName, "in model",
				     omxQuotes(modelname), "does not contain dimnames.")
			stop(msg, call. = FALSE)	
		}
		covRows <- covDimnames[[1]]
		covCols <- covDimnames[[2]]
		if (is.null(covRows) || is.null(covCols) ||
		    (length(covRows) != length(covCols)) || !all(covRows == covCols)) {
			msg <- paste("The",type,"covariance matrix associated",
				     "with", expectationName, "in model",
				     omxQuotes(modelname), "does not contain identical",
				     "row and column dimnames.")
			stop(msg, call.=FALSE)
		}
	}
	if (is.null(means) || (!isS4(means) && is.na(means)) || !length(means)) return()
	meanDimnames <- dimnames(means)
	if (is.null(meanDimnames)) {
			msg <- paste("The",type,"means matrix associated",
				"with", expectationName, "in model",
				omxQuotes(modelname), "does not contain dimnames.")
			stop(msg, call.=FALSE)	
	}
	meanRows <- meanDimnames[[1]]
	meanCols <- meanDimnames[[2]]
	meanNames <- c()
	if (length(meanRows) > length(meanCols)) { meanNames <- meanRows } else { meanNames <- meanCols }
	if ((length(covCols) != length(meanNames)) || !all(covCols == meanNames)) {
			msg <- paste("The",type,"covariance and",type,
				"means matrices associated",
				"with", expectationName, "in model",
				omxQuotes(modelname), "do not contain identical",
				"dimnames.")
			stop(msg, call.=FALSE)
	}
}

verifyObservedNames <- function(data, means, type, flatModel, modelname, expectationName) {
	dataNames <- dimnames(data)
	if(is.null(dataNames)) {
		msg <- paste("The observed data associated with the",
			expectationName, "expectation function in model",
			omxQuotes(modelname), "does not contain dimnames.")
		stop(msg, call. = FALSE)
	}
	if (type == "cov" || type == "cor") {
		if (length(dataNames) < 2 ||
			is.null(dataNames[[1]]) || is.null(dataNames[[2]]) || 
			!identical(dataNames[[1]], dataNames[[2]])) {
				msg <- paste("The dataset associated with the", expectationName,
					"expectation function in model", omxQuotes(modelname),
    	            "does not contain identical row and column non-NULL dimnames.")
			stop(msg, call. = FALSE)
		}
		if (!single.na(means) && is.null(dimnames(means))) {
			msg <- paste("In model", omxQuotes(modelname), 
				", the observed means vector does not contain column names.",
				"Use the name() function to assign names to the means vector.")
			stop(msg, call. = FALSE)
		}
		if (!single.na(means) && !identical(dataNames[[1]], dimnames(means)[[2]])) {
			msg <- paste("The observed covariance or correlation matrix associated with the", expectationName,
				"expectation function in model", omxQuotes(modelname),
				"does not contain identical dimnames to the observed means vector.")
			stop(msg, call. = FALSE)
		}
	} else if ((type == "raw") && (length(dataNames) < 2 || is.null(dataNames[[2]]))) {
		msg <- paste("The dataset associated with the", expectationName,
				"expectation function in model", omxQuotes(modelname),
				"does not contain column names (use dimnames).")
		stop(msg, call. = FALSE)
	}
}

generateDataColumns <- function(flatModel, covNames, dataName) {
	retval <- c()
	if (length(covNames) == 0) return(retval)
	dataColumnNames <- colnames(flatModel@datasets[[dataName]]@observed)
	for(i in 1:length(covNames)) {
		targetName <- covNames[[i]]
		index <- match(targetName, dataColumnNames)
		if(is.na(index)) {
			msg <- paste("The column name", omxQuotes(targetName),
				"in the expected covariance matrix",
				"of the expectation function in model",
				omxQuotes(flatModel@name),
				"cannot be found in the column names of the data.")
			stop(msg, call. = FALSE)
		}
		retval[[i]] <- index - 1
	}
	return(retval)
}


updateThresholdDimnames <- function(flatExpectation, flatModel, labelsData) {
	threshName <- flatExpectation@thresholds
	if (is.na(threshName)) {
		return(flatModel)
	}
	thresholds <- flatModel[[threshName]]
	if (is.null(thresholds)) {
		modelname <- getModelName(flatExpectation)
		stop(paste("Unknown thresholds name", 
			omxQuotes(simplifyName(threshName, modelname)),
			"detected in the expectation function",
			"of model", omxQuotes(modelname)), call. = FALSE)
	}
	dims <- flatExpectation@threshnames
	if (!is.null(colnames(thresholds)) && !single.na(dims) && 
		!identical(colnames(thresholds), dims)) {
		modelname <- getModelName(flatExpectation)
		msg <- paste("The thresholds matrix associated",
		"with the expectation function in model", 
		omxQuotes(modelname), "contains column names and",
		"the expectation function has specified non-identical threshnames.")
		stop(msg, call.=FALSE)      
	}
	if (is.null(colnames(thresholds)) && !single.na(dims)) {
		tuple <- evaluateMxObject(threshName, flatModel, labelsData, new.env(parent = emptyenv()))
		threshMatrix <- tuple[[1]]
		if (ncol(threshMatrix) != length(dims)) {
			modelname <- getModelName(flatExpectation)
			msg <- paste("The thresholds matrix associated",
			"with the expectation function in model", 
			omxQuotes(modelname), "is not of the same length as the 'threshnames'",
			"argument provided by the expectation function. The 'threshnames' argument is",
			"of length", length(dims), "and the expected covariance matrix",
			"has", ncol(threshMatrix), "columns.")
			stop(msg, call.=FALSE)      
		}
		dimnames(flatModel[[threshName]]) <- list(NULL, dims)
	}
	return(flatModel)
}

updateExpectationDimnames <- function(flatExpectation, flatModel,
		labelsData, unsafe = FALSE) {
	covName <- flatExpectation@covariance
	meansName <- flatExpectation@means
	if (is.na(meansName)) {
		means <- NA
	} else {
		means <- flatModel[[meansName]]
	}
	covariance <- flatModel[[covName]]
	if (is.null(covariance)) {
		modelname <- getModelName(flatExpectation)
		stop(paste("Unknown expected covariance name", 
			omxQuotes(simplifyName(covName, modelname)),
			"detected in the expectation function",
			"of model", omxQuotes(modelname)), call. = FALSE)
	}
	if (is.null(means)) {
		modelname <- getModelName(flatExpectation)
		stop(paste("Unknown expected means name", 
			omxQuotes(simplifyName(meansName, modelname)),
			"detected in the expectation function",
			"of model", omxQuotes(modelname)), call. = FALSE)
	}
	dims <- flatExpectation@dims
	if (!is.null(dimnames(covariance)) && !single.na(dims) && 
		!identical(dimnames(covariance), list(dims, dims))) {
		modelname <- getModelName(flatExpectation)
		msg <- paste("The expected covariance matrix associated",
			"with the expectation function in model", 
			omxQuotes(modelname), "contains dimnames: ", 
            paste(toString(dimnames(covariance)), ".", sep = ""),
			"The expectation function has specified dimnames:", 
			paste(toString(dims), ".", sep =""))
		stop(msg, call.=FALSE)		
	}
	if (is.null(dimnames(covariance)) && !single.na(dims)) {
		if (!unsafe) {
			tuple <- evaluateMxObject(covName, flatModel, labelsData, new.env(parent = emptyenv()))
			covMatrix <- tuple[[1]]
			if (nrow(covMatrix) != ncol(covMatrix)) {
				modelname <- getModelName(flatExpectation)
				msg <- paste("The expected covariance matrix associated",
					"with the expectation function in model", 
					omxQuotes(modelname), "is not a square matrix.")
				stop(msg, call.=FALSE)		
			}
			if (nrow(covMatrix) != length(dims)) {
				modelname <- getModelName(flatExpectation)
				msg <- paste("The expected covariance matrix associated",
					"with the expectation function in model", 
					omxQuotes(modelname), "is not of the same length as the 'dimnames'",
					"argument provided by the expectation function. The 'dimnames' argument is",
					"of length", length(dims), "and the expected covariance matrix",
					"has", nrow(covMatrix), "rows and columns.")
				stop(msg, call.=FALSE)		
			}
		}
		dimnames(flatModel[[covName]]) <- list(dims, dims)
	}

	if (!isS4(means) && is.na(means)) {
		return(flatModel)
	}

	if (!is.null(dimnames(means)) && !single.na(dims) &&
		!identical(dimnames(means), list(NULL, dims))) {
		modelname <- getModelName(flatExpectation)
		msg <- paste("The expected means matrix associated",
			"with the expectation function in model", 
			omxQuotes(modelname), "contains dimnames: ", 
            paste(toString(dimnames(means)), ".", sep = ""),
			"The expectation function has specified dimnames:", 
			paste(toString(dims), ".", sep =""))
		stop(msg, call.=FALSE)	
	}
	if (is.null(dimnames(means)) && !single.na(dims)) {
		if (!unsafe) {
			tuple <- evaluateMxObject(meansName, flatModel, labelsData, new.env(parent = emptyenv()))
			meansMatrix <- tuple[[1]]
			if (nrow(meansMatrix) != 1) {
				modelname <- getModelName(flatExpectation)
				msg <- paste("The expected means vector associated",
					"with the expectation function in model", 
					omxQuotes(modelname), "is not a 1 x n matrix.",
					"It has dimensions", nrow(meansMatrix), "x", 
					paste(ncol(meansMatrix), '.', sep=''))
				stop(msg, call.=FALSE)		
			}
			if (ncol(meansMatrix) != length(dims)) {
				modelname <- getModelName(flatExpectation)
				msg <- paste("The expected means vector associated",
					"with the expectation function in model", 
					omxQuotes(modelname), "is not of the same length as the 'dimnames'",
					"argument provided by the expectation function. The 'dimnames' argument is",
					"of length", length(dims), "and the expected means vector",
					"has", ncol(meansMatrix), "columns.")
				stop(msg, call.=FALSE)
			}
		}
		dimnames(flatModel[[meansName]]) <- list(NULL, dims)
	}
	return(flatModel)
}

checkThreshnames <- function(threshnames) {
	if (single.na(threshnames)) threshnames <- as.character(NA)
	if (!is.vector(threshnames) || typeof(threshnames) != 'character') {
		stop("'threshnames' argument is not a character vector")
	}
	if (length(threshnames) == 0) {
		stop("'threshnames' argument cannot be an empty vector")
	}
	if (length(threshnames) > 1 && any(is.na(threshnames))) {
		stop("NA values are not allowed for 'threshnames' vector")
	}
	tt <- table(threshnames)
	if (any(tt > 1)) {
		stop(paste("'threshnames' argument contains", omxQuotes(names(tt)[tt > 1]),
			   "more than once. \nIf you are having problems with Doppelgangers",
			   "perhaps you should check the basement for pods :)"))
	}
	return(threshnames)
}

mxExpectationNormal <- function(covariance, means = NA, 
	dimnames = NA, thresholds = NA, threshnames = dimnames) {
	if (missing(covariance) || typeof(covariance) != "character") {
		stop("'covariance' argument is not a string (the name of the expected covariance matrix)")
	}
	if (!(single.na(means) || typeof(means) == "character")) {
		stop("Means argument is not a string (the name of the expected means matrix)")
	}
	if (is.na(means)) means <- as.integer(NA)
	if (single.na(thresholds)) thresholds <- as.character(NA)
	if (single.na(dimnames)) dimnames <- as.character(NA)
	if (!is.vector(dimnames) || typeof(dimnames) != 'character') {
		stop("'dimnames' argument is not a character vector")
	}
	if (length(thresholds) != 1) {
		stop("'thresholds' argument must be a single matrix or algebra name")
	}
	if (length(dimnames) == 0) {
		stop("'dimnames' argument cannot be an empty vector")
	}
	if (length(dimnames) > 1 && any(is.na(dimnames))) {
		stop("NA values are not allowed for 'dimnames' vector")
	}
	threshnames <- checkThreshnames(threshnames)
	return(new("MxExpectationNormal", covariance, means, dimnames, thresholds, threshnames))
}

displayMxExpectationNormal <- function(expectation) {
	cat("MxExpectationNormal", omxQuotes(expectation@name), '\n')
	cat("$covariance :", omxQuotes(expectation@covariance), '\n')
	cat("$means :", omxQuotes(expectation@means), '\n')
	if (single.na(expectation@dims)) {
		cat("$dims : NA \n")
	} else {
		cat("$dims :", omxQuotes(expectation@dims), '\n')
	}
	if (single.na(expectation@thresholds)) {
		cat("$thresholds : NA \n")
	} else {
		cat("$thresholds :", omxQuotes(expectation@thresholds), '\n')
	}
	if (single.na(expectation@threshnames)) {
		cat("$threshnames : NA \n")
	} else {
		cat("$threshnames :", omxQuotes(expectation@threshnames), '\n')
	}
	invisible(expectation)
}


setMethod("print", "MxExpectationNormal", function(x,...) { 
	displayMxExpectationNormal(x) 
})

setMethod("show", "MxExpectationNormal", function(object) { 
	displayMxExpectationNormal(object) 
})
