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

##' MxRAMModel
##'
##' This is an internal class and should not be used directly.
##'
##' @aliases
##' $<-,MxRAMModel-method
##' [[<-,MxRAMModel-method
setClass(Class = "MxRAMModel",
	representation = representation(),
	contains = "MxModel")

imxModelTypes[['RAM']] <- "MxRAMModel"

##' imxVariableTypes
##'
##' This is an internal function exported for those people who know
##' what they are doing.
##'
##' @details The acceptable variable types
imxVariableTypes <- c(imxVariableTypes, "exogenous", "endogenous")

# Define generic functions

setMethod("imxTypeName", "MxRAMModel", 
	function(model) { "RAM" }
)

setMethod("imxInitModel", "MxRAMModel", 
	function(model) {
		if (any(!is.na(match(model@name, c('A','S','F','M'))))) {
			stop(paste("The name", omxQuotes(model@name), "is not valid",
				   "for a type='RAM' model:\n Pick something other than ",
				   "reserved names ", omxQuotes(c('A','S','F','M'))))
		}
		if (is.null(model[['A']])) {
			model[['A']] <- createMatrixA(model)
		}
		if (is.null(model[['S']])) {
			model[['S']] <- createMatrixS(model)
		}
		if (is.null(model[['expectation']])) {
			model[['expectation']] <- mxExpectationRAM('A', 'S', 'F')
		}
		if (is.null(model[['fitfunction']])) {
			model[['fitfunction']] <- mxFitFunctionML()
		}
		model[['F']] <- createMatrixF(model)
		return(model)
	}
)

setMethod("imxModelBuilder", "MxRAMModel", 
	function(model, lst, name, 
		 manifestVars, latentVars, submodels, remove, independent) {
		model <- nameArgument(model, name)
		model <- variablesArgumentRAM(model, manifestVars, latentVars, submodels, remove)
		model <- listArgumentRAM(model, lst, remove)
		notPathOrData <- getNotPathsOrData(lst)
		callNextMethod(model, notPathOrData, NA, character(), 
			character(), list(), remove, independent)
	}
)

setMethod("imxVerifyModel", "MxRAMModel",
	  function(model) {
		  if ((length(model$A) == 0) ||
		      (length(model$S) == 0)) {
			  msg <- paste("The RAM model", omxQuotes(model@name),
				       "does not contain any paths.",
				       " You can add paths to your model like this:",
				       " mxPath(from = 'x1', to = 'y1')")
			  stop(msg, call. = FALSE)
		  }
		  expectation <- model$expectation
		  if (!is.null(expectation) && is(expectation, "MxExpectationRAM")) {
			  if (!is.null(model@data) && !single.na(model@data@means) &&
			      is.null(model$M)) {
				  msg <- paste("The RAM model", omxQuotes(model@name),
					       "contains an observed means vector",
					       "but has not specified any means paths.")
				  stop(msg, call. = FALSE)
			  }
			  if (!is.null(model@data)) {
				  threshNames <- intersect(getDataThresholdNames(model@data), model@manifestVars)
				  # Only pay attention to (1) manifest variables (2) that need thresholds.
				  # This saves the case where an mxFactor is used as a definition variable.
				  if(!is.null(model[[expectation@thresholds]])) {
					  missingThresholds <- setdiff(threshNames, colnames(model[[expectation@thresholds]]))
				  } else {
					  missingThresholds <- threshNames
				  }
				  if(length(missingThresholds)) {
					  msg <- paste("The RAM model", omxQuotes(model@name),
						       "contains data that requires thresholds for columns",
						       omxQuotes(missingThresholds), "but has not specified any",
							   "thresholds for those columns.",
							   "You can specify thresholds for your model like this:",
							   "mxThreshold(vars='x1', nThresh=1, values=0)")
					  stop(msg, call. = FALSE)
				  }
			  }
		  }
		  if (length(model@submodels) > 0) {
			  return(all(sapply(model@submodels, imxVerifyModel)))
		  }
		  return(TRUE)
	  })


setReplaceMethod("[[", "MxRAMModel",
	function(x, i, j, value) {
		return(replaceMethodRAM(x, i, value))
	}
)

setReplaceMethod("$", "MxRAMModel",
	function(x, name, value) {
		return(replaceMethodRAM(x, name, value))
	}
)

# Helper functions used by the generic functions

variablesArgumentRAM <- function(model, manifestVars, latentVars, submodels, remove) {
	if (single.na(manifestVars)) {
		manifestVars <- character()
	}
	if (single.na(latentVars)) {
		latentVars <- character()
	}
	if (remove == TRUE) {
		if (length(latentVars) + length(manifestVars) > 0) {
			model <- removeVariablesRAM(model, latentVars, manifestVars)
		}
		if (length(submodels)) for(i in 1:length(submodels)) {
			model <- removeSingleNamedEntity(model, submodels[[i]])
		}
	} else {
		if (length(manifestVars) + length(latentVars) > 0) {
			latentVars <- varsToCharacter(latentVars, "latent")
			manifestVars <- varsToCharacter(manifestVars, "manifest")
			checkVariables(model, latentVars, manifestVars)
			model <- addVariablesRAM(model, latentVars, manifestVars)
		}
		if (length(submodels)) for(i in 1:length(submodels)) {
			model <- addSingleNamedEntity(model, submodels[[i]])
		}
	}
	return(model)
}

removeVariablesRAM <- function(model, latent, manifest) {
	missingLatent <- setdiff(latent, model@latentVars)
	missingManifest <- setdiff(manifest, model@manifestVars)
	if (length(missingLatent) > 0) {
		stop(paste("The latent variable(s)", omxQuotes(missingLatent),
			"are not present in the model.",
			"They cannot be deleted"), call. = FALSE)
	} else if (length(missingManifest) > 0) {
		stop(paste("The manifest variable(s)", omxQuotes(missingManifest),
			"are not present in the model.",
			"They cannot be deleted"), call. = FALSE)
	} else if (length(unique(latent)) != length(latent)) {
		stop("The latent variables list contains duplicate elements",
			call. = FALSE)
	} else if (length(unique(manifest)) != length(manifest)) {
		stop("The manifest variables list contains duplicate elements",
			call. = FALSE)
	}
	model@latentVars <- setdiff(model@latentVars, latent)
	model@manifestVars <- setdiff(model@manifestVars, manifest)
	A <- model[['A']]
	S <- model[['S']]
	A <- removeVariablesAS(A, latent)
	A <- removeVariablesAS(A, manifest)
	S <- removeVariablesAS(S, latent)
	S <- removeVariablesAS(S, manifest)
	model[['F']] <- createMatrixF(model)
	model[['A']] <- A
	model[['S']] <- S
	Thresh <- model[['Thresholds']]
	if(!all.na(Thresh)) {
		newCols <- setdiff(colnames(Thresh), manifest)
		newRows <- nrow(Thresh) - min(colSums(is.na(Thresh@values[,newCols])))
		model[['Thresholds']] <- Thresh[1:newRows,newCols]
	}
	return(model)
}

addVariablesRAM <- function(model, latent, manifest) {

	modelLatent   <- unlist(model@latentVars, use.names = FALSE)
	modelManifest <- unlist(model@manifestVars, use.names = FALSE)

	model <- addVariablesHelper(model, "latentVars", latent)
	model <- addVariablesHelper(model, "manifestVars", manifest)

	latent <- unlist(latent, use.names = FALSE)
	manifest <- unlist(manifest, use.names = FALSE)

	newLatent   <- setdiff(latent, modelLatent)
	newManifest <- setdiff(manifest, modelManifest)

	A <- model[['A']]
	S <- model[['S']]
	M <- model[['M']]
	if (is.null(A)) {
		A <- createMatrixA(model)
	} else {
		A <- addVariablesAS(A, model, newLatent, newManifest)
	}
	if (is.null(S)) {
		S <- createMatrixS(model)
	} else {
		S <- addVariablesAS(S, model, newLatent, newManifest)
	}
	if (!is.null(M)) {
		M <- addVariablesM(M, model, newLatent, newManifest)
	}
	model[['A']] <- A
	model[['S']] <- S
	model[['M']] <- M
	model[['F']] <- createMatrixF(model)
	return(model)
}


listArgumentRAM <- function(model, lst, remove) {
	if(remove == TRUE) {
		model <- removeEntriesRAM(model, lst)
	} else {
		model <- addEntriesRAM(model, lst)
	}
	return(model)
}

addEntriesRAM <- function(model, entries) {
	if (length(entries) == 0) {
		return(model)
	}

	filter <- sapply(entries, is, "MxPath")
	paths <- entries[filter]
	if (length(paths) > 0) {
		model <- insertAllPathsRAM(model, paths)
	}
	filter <- sapply(entries, is, "MxThreshold")
	thresholds <- entries[filter]
	if(length(thresholds) > 0) {
		model <- insertAllThresholdsRAM(model, thresholds)
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
		model[['F']] <- createMatrixF(model)
	}
	return(model)
}

requireThresholds <- function(data) {
	return(!is.null(data) && ((data@type == 'raw') ||
		((data@type == 'cov' || data@type == 'cor'))))
}

requireMeansVector <- function(data) {
	return(!is.null(data) && ((data@type == 'raw') ||
		((data@type == 'cov' || data@type == 'cor') &&
		 !(length(data@means) == 1 && is.na(data@means)))))
}

removeEntriesRAM <- function(model, entries) {
	if (length(entries) == 0) {
		return(model)
	}
	filter <- sapply(entries, is, "MxPath")
	paths <- entries[filter]
	if (length(paths) > 0) {
		model <- removeAllPathsRAM(model, paths)
	}
	return(model)
}

getNotPathsOrData <- function(lst) {
	if (length(lst) == 0) {
		return(lst)
	}
	pathfilter <- sapply(lst, is, "MxPath")
	datafilter <- sapply(lst, is, "MxData")
	thresholdfilter <- sapply(lst, is, "MxThreshold")
	retval <- lst[!(pathfilter | datafilter | thresholdfilter)]
	return(retval)
}

expectationIsMissingThresholds <- function(model) {
	expectation <- model@expectation
	return(!is.null(expectation) &&
	is(expectation, "MxExpectationRAM") &&
		is.na(expectation@thresholds))
}


expectationIsMissingMeans <- function(model) {
	expectation <- model@expectation
	return(!is.null(expectation) &&
		is(expectation, "MxExpectationRAM") &&
		is.na(expectation@M))
}

insertAllThresholdsRAM <- function(model, thresholds) {
	Thresh <- model[[model@expectation@thresholds]]
	if (is.null(Thresh)) { 
		Thresh <- mxMatrix("Full", 0, 0, name="Thresholds", condenseSlots=FALSE)
		if(expectationIsMissingThresholds(model)) {
			model@expectation@thresholds <- "Thresholds"
		} else {
			Thresh <- model[[model@expectation@thresholds]]
		}
	}
	
	legalVars <- model@manifestVars
	isUsed <- matrix(FALSE, 1, length(legalVars))
	colnames(isUsed) <- legalVars
	isUsed[colnames(Thresh)] <- TRUE
	maxNThresh <- nrow(Thresh)

	allVars <- unique(as.character(lapply(thresholds, getElement, "variable")))
	varExist <- allVars %in% legalVars
	if(!all(varExist)) {
		missingVars <- allVars[!varExist]
		stop(paste("Nice try, you need to add", 
		omxQuotes(missingVars), 
			"to the manifestVars before you",
			"can assign them thresholds."), call. = FALSE)
	}

	maxNThresh <- max(sapply(thresholds, getElement, "nThresh"))
	isUsed[allVars] <- TRUE
	
	newVars <- union(colnames(Thresh), allVars)
	if(length(newVars) > ncol(Thresh)) {  # Rebuild Threshold matrix if needed
		oldCols <- ncol(Thresh)
		oldRows <- nrow(Thresh)
		newCols <- length(newVars)
		newRows <- max(nrow(Thresh), maxNThresh)
		newThresh <- mxMatrix("Full", newRows, newCols, dimnames=list(NULL,
			newVars), name=Thresh@name, condenseSlots=FALSE)  # Maintains the old ordering
		if(oldRows > 0 && oldCols > 0) {
			newThresh[1:oldRows, 1:oldCols] <- Thresh
		}
		Thresh <- newThresh 
	}
	if(!is.list(thresholds)) { thresholds <- list(thresholds)}
	for(i in 1:length(thresholds)) {
		thisThresh <- thresholds[[i]]
		values <- thisThresh@values
		msg = NULL
		if(all.na(values)) {
			msg = paste("The thresholds you are attempting to specify",
				"does not have any starting values, but type='RAM' models require them.")
		}
		if(!identical(values, sort(values, na.last=NA))) { 
			msg = paste("The thresholds you are attempting to specify",
				"has starting values that are not strictly increasing,",
				"but type='RAM' models require them to be.")
		}
		if(length(msg)>0) {
			msg <- paste(msg, paste("An easy way to specify threshold starting values",
				"that are evenly spaced across a normal distribution is using:"),
				"mxThreshold(vars='x1', nThresh=3, values=mxNormalQuantiles(3))",
				"See '?mxNormalQuantiles' and '?mxThreshold' for more details.", sep='\n')
				stop(msg, .call=FALSE)
		}
		thisVar <- thisThresh@variable
		theseRows <- 1:thisThresh@nThresh
		threshMat <- as(thisThresh, "MxMatrix")
		Thresh[theseRows, thisVar] <- threshMat
	}

	model[['Thresholds']] <- Thresh
	
	return(model)
}

insertAllPathsRAM <- function(model, paths) {
	A <- model[['A']]
	S <- model[['S']]
	M <- model[['M']]
	if (is.null(A)) { A <- createMatrixA(model) }
	if (is.null(S)) { S <- createMatrixS(model) }
	
	legalVars <- c(model@latentVars, model@manifestVars, "one")
	
	for(i in 1:length(paths)) {
		path <- paths[[i]]
	
		missingvalues <- is.na(path@values)
		path@values[missingvalues] <- 0
		
		if (single.na(path@to)) {
			path@to <- path@from
			paths[[i]] <- path
		}
		
		if (!is.na(path@joinKey)) {
			upperFrom <- strsplit(path@from, imxSeparatorChar, fixed = TRUE)
			upperFromBadLen <- sapply(upperFrom, length) != 2
			if (any(upperFromBadLen)) {
				msg <- paste("Between level paths must be from",
					     "modelName.variableName, not",
					     omxQuotes(path@from[upperFromBadLen]))
				stop(msg, call.=FALSE)
			}
			fromUpperModels <- sapply(upperFrom, function(v) v[1])
			if (length(unique(fromUpperModels)) > 1) {
				msg <- paste("Deal with one upper level model at a time, not",
					     omxQuotes(unique(fromUpperModels)))
				stop(msg, call.=FALSE)
			}
			fromUpperVars <- sapply(upperFrom, function(v) v[2])

			upperModelName <- fromUpperModels[1]
			upperModel <- model[[upperModelName]]
			if (is.null(upperModel)) {
				msg <- paste("Nice try. You need to create an upper level RAM",
				      "model called", omxQuotes(upperModelName),
				      "and add it as a submodel of", omxQuotes(model@name),
				      "before you can create paths between these models.")
				stop(msg, call.=FALSE)
			}

			upperVars <- c(upperModel@manifestVars, upperModel@latentVars)
			upperVarExist <- fromUpperVars %in% upperVars
			if (!all(upperVarExist)) {
				stop(paste("Nice try, you need to add",
					   omxQuotes(fromUpperVars[!upperVarExist]),
					   "to either manifestVars or latentVars in model",
					   omxQuotes(upperModelName),
					   "before you can use them in a path."), call. = FALSE)
			}

			lowerVarExist <- path@to %in% rownames(A)
			if (!all(lowerVarExist)) {
				stop(paste("Nice try, you need to add",
					   omxQuotes(path@to[!lowerVarExist]),
					   "to either manifestVars or latentVars before you",
					   "can use them in a path."), call. = FALSE)
			}

			bMatName <- NULL
			priorBetween <- model$expectation$between
			sameJoinMask <- sapply(priorBetween, function(x) {
				bmat <- model[[ x ]]
				if (is.null(bmat)) return(FALSE)
				bmat$joinKey == path@joinKey && bmat$joinModel == upperModelName
			})
			if (length(sameJoinMask) && sum(sameJoinMask) > 1) {
				stop(paste("Confusingly there is more than 1 join to", omxQuotes(upperModelName),
					   "using foreign key", omxQuotes(path@joinKey)), call.=FALSE)
			}
			if (all(!sameJoinMask)) {
				bMatNameBase <- paste0("from_", upperModelName)
				bMatName <- bMatNameBase
				found <- FALSE
				for (try in 1:9) {
					if (is.null( model[[ bMatName ]] )) {
						found = TRUE
						break
					}
					bMatName <- paste0(bMatNameBase, try)
				}
				if (!found) {
					stop(paste("Failed to invent an unused name for the",
						   "between level mapping matrix. Tried variations on",
						   omxQuotes(bMatNameBase)), call.=FALSE)
				}
				model <- mxModel(model, mxMatrix(name=bMatName, nrow=nrow(A), ncol=length(upperVars),
								 dimnames=list(rownames(A), upperVars),
								 joinKey = path@joinKey, joinModel = upperModelName))
			} else {
				bMatName <- priorBetween[ which(sameJoinMask) ]
			}

			bMat <- model[[ bMatName ]]
			bMat <- imxGentleResize(bMat, list(rownames(A), upperVars))

			expanded <- expandPathConnect(fromUpperVars, path@to, path@connect)
			allfrom <- expanded$from
			allto   <- expanded$to

			maxlength <- max(length(path@from), length(path@to)) - 1L
			for(i in 0:maxlength) {
				from <- allfrom[[i %% length(allfrom) + 1L]]
				to <- allto[[i %% length(allto) + 1L]]
				bMat$values[to, from] <- path@values[[i %% length(path@values) + 1L]]
				bMat$free  [to, from] <- path@free  [[i %% length(path@free)   + 1L]]
				bMat$labels[to, from] <- path@labels[[i %% length(path@labels) + 1L]]
				bMat$ubound[to, from] <- path@ubound[[i %% length(path@ubound) + 1L]]
				bMat$lbound[to, from] <- path@lbound[[i %% length(path@lbound) + 1L]]
			}

			model[[ bMatName ]] <- bMat
			model$expectation$between <- unique(c(priorBetween, bMatName))
			next
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
			if (is.null(M)) {
				M <- createMatrixM(model) 
				if(expectationIsMissingMeans(model)) {
					model@expectation@M <- "M"
				}
			}
			M <- insertMeansPathRAM(path, M)
		} else {
			if (length(intersect(path@from, "one"))) {
				# This really ought to work, but it doesn't.
				msg <- paste("Cannot create path from 'one' and other nodes simultaneously.",
					     "Create paths from 'one' and then separately from",
					     omxQuotes(setdiff(path@from, "one")))
				stop(msg, call.=FALSE)
			}
			expanded <- expandPathConnect(path@from, path@to, path@connect)
			path@from <- expanded$from
			path@to   <- expanded$to
			retval <- insertPathRAM(path, A, S)
			A <- retval[[1]]
			S <- retval[[2]]	
		}
	}
	model[['A']] <- A
	model[['S']] <- S
	if (!is.null(M)) {
		model[['M']] <- M
	}
	
	return(model)
}

removeAllPathsRAM <- function(model, paths) {
	A <- model[['A']]
	S <- model[['S']]
	M <- model[['M']]
	if (is.null(A)) { A <- createMatrixA(model) }
	if (is.null(S)) { S <- createMatrixS(model) }
	for(i in 1:length(paths)) {

		path <- paths[[i]]

		if (single.na(path@to)) {
			path@to <- path@from
			paths[[i]] <- path
		}
		
		if (length(path@from) == 1 && (path@from == "one")) {		
			M <- removeMeansPathRAM(path, M)
		} else {
			expanded <- expandPathConnect(path@from, path@to, path@connect)
			path@from <- expanded$from
			path@to   <- expanded$to
			retval <- removePathRAM(path, A, S)
			A <- retval[[1]]
			S <- retval[[2]]
		}
	}
	model[['A']] <- A
	model[['S']] <- S
	if (!is.null(M)) {
		model[['M']] <- M
	}
	return(model)
}


insertPathRAM <- function(path, A, S) {
	allfrom <- path@from
	allto <- path@to
	allarrows <- path@arrows
	allfree <- path@free
	allvalues <- path@values
	alllabels <- path@labels
	alllbound <- path@lbound
	allubound <- path@ubound
	maxlength <- max(length(allfrom), length(allto))
	A_free <- A@free
	A_values <- A@values
	A_labels <- A@labels
	A_lbound <- A@lbound
	A_ubound <- A@ubound
	S_free   <- S@free
	S_values <- S@values
	S_labels <- S@labels
	S_lbound <- S@lbound
	S_ubound <- S@ubound
	for(i in 0:(maxlength - 1)) {
		from <- allfrom[[i %% length(allfrom) + 1]]
		to <- allto[[i %% length(allto) + 1]]
		arrows <- allarrows[[i %% length(allarrows) + 1]]
		nextvalue <- allvalues[[i %% length(allvalues) + 1]]
		nextfree <- allfree[[i %% length(allfree) + 1]]
		nextlabel <- alllabels[[i %% length(alllabels) + 1]]
		nextubound <- allubound[[i %% length(allubound) + 1]]
		nextlbound <- alllbound[[i %% length(alllbound) + 1]]		
		if (arrows == 1) {
			A_free[to, from] <- nextfree
			A_values[to, from] <- nextvalue
			A_labels[to, from] <- nextlabel
			A_ubound[to, from] <- nextubound
			A_lbound[to, from] <- nextlbound
			S_values[to, from] <- 0
			S_labels[to, from] <- as.character(NA)
			S_free[to, from] <- FALSE
			S_values[from, to] <- 0
			S_labels[from, to] <- as.character(NA)
			S_free[from, to] <- FALSE			
		} else if (arrows == 2) {
			S_free[to, from] <- nextfree
			S_values[to, from] <- nextvalue
			S_labels[to, from] <- nextlabel
			S_ubound[to, from] <- nextubound
			S_lbound[to, from] <- nextlbound
			S_free[from, to] <- nextfree
			S_values[from, to] <- nextvalue
			S_labels[from, to] <- nextlabel
			S_ubound[from, to] <- nextubound
			S_lbound[from, to] <- nextlbound
			A_values[to, from] <- 0
			A_labels[to, from] <- as.character(NA)
			A_free[to, from] <- FALSE
			A_values[to, from] <- 0
			A_labels[to, from] <- as.character(NA)
			A_free[to, from] <- FALSE
		} else {
			stop(paste("Unknown arrow type", arrows, 
					"with source", omxQuotes(from), 
				"and sink", omxQuotes(to)),
				call. = FALSE)
		}
	}
	A@free <- A_free
	A@values <-	A_values 
	A@labels <-	A_labels 
	A@lbound <-	A_lbound 
	A@ubound <-	A_ubound 
	S@free <- S_free   
	S@values <-	S_values 
	S@labels <-	S_labels 
	S@lbound <- S_lbound 
	S@ubound <- S_ubound 
	return(list(A, S))
}

removeMeansPathRAM <- function(path, M) {
	if(is.null(M)) {
		return(NULL)
	}
	allto <- path@to
	for(i in 0:(length(allto) - 1)) {
		to <- allto[[i %% length(allto) + 1]]
		M@free[1, to] <- FALSE
		M@values[1, to] <- 0
		M@labels[1, to] <- as.character(NA)
	}
	return(M)
}

insertMeansPathRAM <- function(path, M) {
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
		M@free[1, to] <- nextfree
		M@values[1, to] <- nextvalue
		M@labels[1, to] <- nextlabel
		M@ubound[1, to] <- nextubound
		M@lbound[1, to] <- nextlbound
	}
	return(M)
}

removePathRAM <- function(path, A, S) {
	allfrom <- path@from
	allto <- path@to
	allarrows <- path@arrows
	maxlength <- max(length(allfrom), length(allto))
	A_free <- A@free
	A_values <- A@values
	A_labels <- A@labels
	S_free   <- S@free
	S_values <- S@values
	S_labels <- S@labels
	for(i in 0:(maxlength - 1)) {
		from <- allfrom[i %% length(allfrom) + 1]
		to <- allto[i %% length(allto) + 1]
		arrows <- allarrows[i %% length(allarrows) + 1]
		if (arrows == 1) {
			A_values[to, from] <- 0
			A_labels[to, from] <- as.character(NA)
			A_free[to, from] <- FALSE		
		} else if (arrows == 2) {
			S_values[to, from] <- 0
			S_labels[to, from] <- as.character(NA)
			S_free[to, from] <- FALSE		
			S_values[from, to] <- 0
			S_labels[from, to] <- as.character(NA)
			S_free[from, to] <- FALSE					
		} else {
			stop(paste("Unknown arrow type", arrows, 
					"with source", omxQuotes(from), 
					"and sink", omxQuotes(to)),
					call. = FALSE)
		}
	}
	A@free <- A_free
	A@values <-	A_values 
	A@labels <-	A_labels 
	S@free <- S_free   
	S@values <-	S_values 
	S@labels <-	S_labels 	
	return(list(A, S))
}

createMatrixM <- function(model) {
	variables <- c(model@manifestVars, model@latentVars)
	len <- length(variables)
	names <- list(NULL, variables)
	values <- matrix(0, 1, len)
	labels <- matrix(as.character(NA), 1, len)
	free <- matrix(c(rep.int(FALSE, length(model@manifestVars)),
		rep.int(FALSE, length(model@latentVars))), 1, len)
	retval <- mxMatrix("Full", values = values, free = free, labels = labels, name = "M", condenseSlots=FALSE)
	dimnames(retval) <- names
	return(retval)
}

createMatrixA <- function(model) {
	variables <- c(model@manifestVars, model@latentVars)
	len <- length(variables)
	names <- list(variables, variables)
	values <- matrix(0, len, len)
	free <- matrix(FALSE, len, len)
	labels <- matrix(as.character(NA), len, len)
	retval <- mxMatrix("Full", values = values, free = free, labels = labels, name = "A", condenseSlots=FALSE)
	dimnames(retval) <- names
	return(retval)
}

createMatrixS <- function(model) {
	variables <- c(model@manifestVars, model@latentVars)
	len <- length(variables)
	names <- list(variables, variables)
	values <- matrix(0, len, len)
	free <- matrix(FALSE, len, len)
	labels <- matrix(as.character(NA), len, len)
	retval <- mxMatrix("Symm", values = values, free = free, labels = labels, name = "S", condenseSlots=FALSE)
	dimnames(retval) <- names
	return(retval)
}

createMatrixF <- function(model) {
	variables <- c(model@manifestVars, model@latentVars)
	len <- length(variables)
	values <- diag(nrow = length(model@manifestVars), ncol = len)
	names <- list(model@manifestVars, variables)
	if (!is.null(model@data) && (model@data@type != 'raw')) {
		manifestNames <- rownames(model@data@observed)
		extraData <- setdiff(manifestNames, model@manifestVars)
		extraVars <- setdiff(model@manifestVars, manifestNames)
		if (length(extraData) > 0) {
			msg <- paste("The observed data contains the variables:",
				omxQuotes(extraData), "that have not been declared in the",
				"manifest variables.")
			stop(msg, call. = FALSE)
		}
		if (length(extraVars) > 0) {
			msg <- paste("The manifest variables include",
				omxQuotes(extraVars), "that have not been found in the",
				"observed data.")
			stop(msg, call. = FALSE)
		}
		dimnames(values) <- names
		values <- values[manifestNames,]
		names <- list(manifestNames, variables)
	}
	free <- matrix(FALSE, length(model@manifestVars), len)
	labels <- matrix(as.character(NA), length(model@manifestVars), len)
	retval <- mxMatrix("Full", values = values, free = free, labels = labels, name = "F", condenseSlots=TRUE)
	dimnames(retval) <- names
	return(retval)
}

addVariablesAS <- function(oldmatrix, model, newLatent, newManifest) {
	newLatent <- length(newLatent)
	newManifest <- length(newManifest)
	oldmatrix@values <- addVariablesMatrix(oldmatrix@values, 0, 
		model, newLatent, newManifest)
	oldmatrix@free <- addVariablesMatrix(oldmatrix@free, FALSE, 
		model, newLatent, newManifest)
	oldmatrix@labels <- addVariablesMatrix(oldmatrix@labels, as.character(NA), 
		model, newLatent, newManifest)
	oldmatrix@lbound <- addVariablesMatrix(oldmatrix@lbound, as.numeric(NA), 
		model, newLatent, newManifest)
	oldmatrix@ubound <- addVariablesMatrix(oldmatrix@ubound, as.numeric(NA), 
		model, newLatent, newManifest)		
	variables <- c(model@manifestVars, model@latentVars)
	dimnames(oldmatrix) <- list(variables, variables)
	return(oldmatrix)
}

addVariablesM <- function(oldmatrix, model, newLatent, newManifest) {
	oldmatrix@values <- addVariablesMatrixM(oldmatrix@values, 0, 0, model, newLatent, newManifest)
	oldmatrix@free   <- addVariablesMatrixM(oldmatrix@free, FALSE, TRUE, model, newLatent, newManifest)
	oldmatrix@labels <- addVariablesMatrixM(oldmatrix@labels, as.character(NA), as.character(NA),
		model, newLatent, newManifest) 
	oldmatrix@lbound <- addVariablesMatrixM(oldmatrix@lbound, as.numeric(NA), as.numeric(NA), 
		model, newLatent, newManifest)
	oldmatrix@ubound <- addVariablesMatrixM(oldmatrix@ubound, as.numeric(NA), as.numeric(NA), 
		model, newLatent, newManifest)
	dimnames(oldmatrix) <- list(NULL, c(model@manifestVars, model@latentVars))
	return(oldmatrix)
}

removeVariablesAS <- function(oldmatrix, variables) {
	if (length(variables) > 0) {
		for (i in 1:length(variables)) {
			index <- match(variables[[i]], dimnames(oldmatrix)[[1]])
			oldmatrix@values <- oldmatrix@values[-index, -index]
			oldmatrix@free <- oldmatrix@free[-index, -index]
			oldmatrix@labels <- oldmatrix@labels[-index, -index]
			oldmatrix@lbound <- oldmatrix@lbound[-index, -index]
			oldmatrix@ubound <- oldmatrix@ubound[-index, -index]
		}
	}
	return(oldmatrix)
}


addVariablesMatrix <- function(oldmatrix, value, model, newLatent, newManifest) {
	currentManifest <- length(model@manifestVars) - newManifest
	currentLatent <- length(model@latentVars) - newLatent
	newSize <- length(model@manifestVars) + length(model@latentVars)
	if (currentManifest > 0) {
		manifestXmanifest <- oldmatrix[1 : currentManifest, 1 : currentManifest]
	} else {
		manifestXmanifest <- matrix(value, 0, 0)
	}
	if (currentLatent > 0) {
		latentStart <- currentManifest + 1
		latentEnd <- currentManifest + currentLatent
		manifestXlatent <- oldmatrix[1 : currentManifest, latentStart : latentEnd]
		latentXmanifest <- oldmatrix[latentStart : latentEnd, 1 : currentManifest]
		latentXlatent <- oldmatrix[latentStart : latentEnd, latentStart : latentEnd]
	} else {
		manifestXlatent <- matrix(value, 0, 0)
		latentXmanifest <- matrix(value, 0, 0)
		latentXlatent <- matrix(value, 0, 0)
	}
	newtop <- cbind(manifestXmanifest, matrix(value, currentManifest, newManifest),
					manifestXlatent, matrix(value, currentManifest, newLatent)) 
	newtop <- rbind(newtop, matrix(value, newManifest, newSize))
	newbottom <- cbind(latentXmanifest, matrix(value, currentLatent, newManifest),
					latentXlatent, matrix(value, currentLatent, newLatent)) 
	newbottom <- rbind(newbottom, matrix(value, newLatent, newSize))
	newmatrix <- rbind(newtop, newbottom)
	return(newmatrix)
}

addVariablesMatrixM <- function(oldmatrix, newLatentValue, newManifestValue, model, newLatent, newManifest) {
	newManifest <- length(newManifest)
	newLatent <- length(newLatent)
	currentManifest <- length(model@manifestVars) - newManifest
	currentLatent <- length(model@latentVars) - newLatent
	values <- c(oldmatrix[1, 1:currentManifest], 
		rep.int(newManifestValue, newManifest),
		oldmatrix[1, (currentManifest + 1) : (currentLatent + currentManifest)],
		rep.int(newLatentValue, newLatent))
	newmatrix <- matrix(values, 1, length(model@manifestVars) + length(model@latentVars))
	return(newmatrix)
}

replaceMethodRAM <- function(model, index, value) {
	pair <- imxReverseIdentifier(model, index)
	namespace <- pair[[1]]
	name <- pair[[2]]
	if (namespace == model@name && name == "data") {
		model@data <- value
		# model@expectation is likely MxExpectationRAM, but it could
		# have been replaced with something else.
		if (requireMeansVector(value) && is(model@expectation, "MxExpectationRAM")) {
			model@expectation@M <- "M"
		}
	} else {
		model <- imxReplaceMethod(model, index, value)
	}
	return(model)
}

