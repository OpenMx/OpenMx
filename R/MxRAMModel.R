#
#   Copyright 2007-2021 by the individuals mentioned in the source code history
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
	contains = "MxModel")

imxModelTypes[['RAM']] <- "MxRAMModel"

##' imxVariableTypes
##'
##' This is an internal function exported for those people who know
##' what they are doing.
##'
##' @details The acceptable variable types
imxVariableTypes <- c("exogenous", "endogenous")

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
		 manifestVars, latentVars, productVars, submodels, remove, independent) {
		model <- nameArgument(model, name)
		model <- variablesArgumentRAM(model, manifestVars, latentVars, productVars, submodels, remove)
		model <- listArgumentRAM(model, lst, remove)
		notPathOrData <- getNotPathsOrData(lst)
		callNextMethod(model, notPathOrData, NA, character(),
			character(), character(), list(), remove, independent)
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
		  if (length(model@submodels) > 0) {
			  return(all(sapply(model@submodels, imxVerifyModel)))
		  }
		  return(TRUE)
	  })


# Helper functions used by the generic functions

variablesArgumentRAM <- function(model, manifestVars, latentVars, productVars, submodels, remove) {
	manifestVars <- unlist(manifestVars)
	latentVars <- unlist(latentVars)
	if (single.na(manifestVars)) {
		manifestVars <- character()
	}
	if (single.na(latentVars)) {
		latentVars <- character()
	}
	if (single.na(productVars)) {
		productVars <- character()
	}
	latentVars <- c(latentVars, productVars)
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
			model[['expectation']]$isProductNode <- colnames(model$A) %in% productVars
      if (length(productVars) && !is.null(model[['M']])) {
        model[['M']]$values[1,productVars] <- 1
      }
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
	if(!isAllNa(Thresh)) {
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
	filter <- sapply(entries, is, "DiscreteBase")
	discrete <- entries[filter]
	if(length(discrete)) {
		model <- insertDiscreteRAM(model, discrete)
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
		if(is(model@data,"MxDataLegacyWLS") && class(model@fitfunction) %in% "MxFitFunctionML"){
			model[['fitfunction']] <- mxFitFunctionWLS()
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
	discretefilter <- sapply(lst, is, "DiscreteBase")
	retval <- lst[!(pathfilter | datafilter | thresholdfilter | discretefilter)]
	return(retval)
}

expectationIsMissing <- function(model, what) {
	expectation <- model@expectation
	return(!is.null(expectation) &&
           is(expectation, "BaseExpectationNormal") &&
           is.na(slot(expectation, what)))
}

expectationIsMissingMeans <- function(model) {
	expectation <- model@expectation
	return(!is.null(expectation) &&
		is(expectation, "MxExpectationRAM") &&
		is.na(expectation@M))
}

insertDiscreteRAM <- function(model, discrete) {
  DiscSpec <- model@expectation@discreteSpec
	Disc <- model[[model@expectation@discrete]]
	if (is.null(Disc)) {
		Disc <- mxMatrix("Full", 0, 0, name="Discrete", condenseSlots=FALSE)
		if(expectationIsMissing(model, 'discrete')) {
			model@expectation@discrete <- "Discrete"
		} else {
			Disc <- model[[model@expectation@discrete]]
		}
	}

	legalVars <- model@manifestVars
	if (is.list(legalVars)) {
	  # for LISREL
	  legalVars <- legalVars$endogenous
	}

	allVars <- unique(sapply(discrete, slot, "variable"))
	varExist <- allVars %in% legalVars
	if(!all(varExist)) {
		missingVars <- allVars[!varExist]
		stop(paste("You need to add", omxQuotes(missingVars),
               "to the manifestVars before you",
               "use them as discrete indicators"), call. = FALSE)
	}

	if(!is.list(discrete)) discrete <- list(discrete)
  todoSpec <- lapply(discrete, function(x) getSpec(x))
  todo <- lapply(discrete, function(x) as(x, 'MxMatrix'))
  newRows <- max(sapply(todo, nrow))

	newVars <- union(colnames(Disc), allVars)
	if(length(newVars) > ncol(Disc)) {
		newDisc <- mxMatrix("Full", newRows, length(newVars),
                        dimnames=list(NULL, newVars), name=Disc@name, condenseSlots=FALSE)
    newDisc$values <- NA
		if (ncol(Disc)) newDisc[1:nrow(Disc), 1:ncol(Disc)] <- Disc
		Disc <- newDisc
    newSpec <- matrix(NA, 2, length(newVars))
    colnames(newSpec) <- newVars
    newSpec[,colnames(DiscSpec)] <- DiscSpec
    DiscSpec <- newSpec
	}

  for (c1 in 1:length(allVars)) {
    n1 <- allVars[c1]
    DiscSpec[,n1] <- todoSpec[[c1]]
    mat <- todo[[c1]]
    Disc[1:nrow(mat),n1] <- mat
  }

  model@expectation@discreteSpec <- DiscSpec
	model[[Disc@name]] <- Disc
  model
}

insertAllThresholdsRAM <- function(model, thresholds) {
	Thresh <- model[[model@expectation@thresholds]]
	if (is.null(Thresh)) {
		Thresh <- mxMatrix("Full", 0, 0, name="Thresholds", condenseSlots=FALSE)
		if(expectationIsMissing(model, 'thresholds')) {
			model@expectation@thresholds <- "Thresholds"
		} else {
			Thresh <- model[[model@expectation@thresholds]]
		}
	}

	legalVars <- model@manifestVars
	if (is.list(legalVars)) {
	  # for LISREL
	  legalVars <- legalVars$endogenous
	}
	maxNThresh <- nrow(Thresh)

	allVars <- unique(as.character(lapply(thresholds, getElement, "variable")))
	varExist <- allVars %in% legalVars
	if(!all(varExist)) {
		missingVars <- allVars[!varExist]
		stop(paste("You need to add",
		omxQuotes(missingVars),
			"to the manifestVars before you",
			"can assign them thresholds."), call. = FALSE)
	}

	maxNThresh <- max(sapply(thresholds, getElement, "nThresh"))

	newVars <- union(colnames(Thresh), allVars)
	if(length(newVars) > ncol(Thresh)) {  # Rebuild Threshold matrix if needed
		oldCols <- ncol(Thresh)
		oldRows <- nrow(Thresh)
		newCols <- length(newVars)
		newRows <- max(nrow(Thresh), maxNThresh)
		newThresh <- mxMatrix("Full", newRows, newCols, dimnames=list(NULL, newVars),
                          name=Thresh@name, condenseSlots=FALSE)  # Maintains the old ordering
    newThresh@values[,] <- NA
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
		if(isAllNa(values)) {
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
				stop(msg, call.=FALSE)
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
	if (is.null(model[['A']])) { model[['A']] <- createMatrixA(model) }
	if (is.null(model[['S']])) { model[['S']] <- createMatrixS(model) }

	for(i in 1:length(paths)) {
		path <- paths[[i]]

		if ("one" %in% path@from && is.null(model[['M']])) {
			model[['M']] <- createMatrixM(model)
			if(expectationIsMissingMeans(model)) {
				model@expectation@M <- "M"
			}
		}

		model <- insertPathRAM(prepPath(path), model)
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


insertPathRAM <- function(path, model) {
	allfrom <- path@from
	allto <- path@to
	allarrows <- path@arrows
	allfree <- path@free
	allvalues <- path@values
	alllabels <- path@labels
	alllbound <- path@lbound
	allubound <- path@ubound
	alljoinKey <- path@joinKey
	maxlength <- max(length(allfrom), length(allto))
	A <- model[['A']]
	S <- model[['S']]
	M <- model[['M']]
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
  selVec <- model[['selectionVector']]
  selPlan <- model$expectation$selectionPlan

	legalVars <- c(colnames(A), "one")
	isProductNode <- model$expectation$isProductNode
	names(isProductNode) <- colnames(A)

	for(i in 0:(maxlength - 1)) {
		from <- allfrom[[i %% length(allfrom) + 1]]
		to <- allto[[i %% length(allto) + 1]]
		arrows <- allarrows[[i %% length(allarrows) + 1]]
		nextvalue <- allvalues[[i %% length(allvalues) + 1]]
		nextfree <- allfree[[i %% length(allfree) + 1]]
		nextlabel <- alllabels[[i %% length(alllabels) + 1]]
		nextubound <- allubound[[i %% length(allubound) + 1]]
		nextlbound <- alllbound[[i %% length(alllbound) + 1]]
		nextjoinKey <- alljoinKey[[i %% length(alljoinKey) + 1]]
    nextstep <- path@step[[i %% length(path@step) + 1]]

		if (!is.na(nextjoinKey)) {
			upperFrom <- strsplit(from, imxSeparatorChar, fixed = TRUE)
			upperFromBadLen <- sapply(upperFrom, length) != 2
			if (any(upperFromBadLen)) {
				msg <- paste("Between level paths must be from",
					     "modelName.variableName, not",
					     omxQuotes(from[upperFromBadLen]))
				stop(msg, call.=FALSE)
			}
			fromUpperModel <- upperFrom[[1]][1]
			fromUpperVar <- upperFrom[[1]][2]

			upperModelName <- fromUpperModel[1]
			upperModel <- model[[upperModelName]]
			if (is.null(upperModel)) {
				msg <- paste("You need to create an upper level RAM",
				      "model called", omxQuotes(upperModelName),
				      "and add it as a submodel of", omxQuotes(model@name),
				      "before you can create paths between these models.")
				stop(msg, call.=FALSE)
			}

			upperVars <- c(upperModel@manifestVars, upperModel@latentVars)
			upperVarExist <- fromUpperVar %in% upperVars
			if (!all(upperVarExist)) {
				stop(paste("You need to add",
					   omxQuotes(fromUpperVar[!upperVarExist]),
					   "to either manifestVars or latentVars in model",
					   omxQuotes(upperModelName),
					   "before you can use them in a path."), call. = FALSE)
			}

			lowerVarExist <- to %in% rownames(A)
			if (!all(lowerVarExist)) {
				stop(paste("You need to add",
					   omxQuotes(to[!lowerVarExist]),
					   "to either manifestVars or latentVars before you",
					   "can use them in a path."), call. = FALSE)
			}

			bMatName <- NULL
			priorBetween <- model$expectation$between
			sameJoinMask <- sapply(priorBetween, function(x) {
				bmat <- model[[ x ]]
				if (is.null(bmat)) return(FALSE)
				bmat$joinKey == nextjoinKey && bmat$joinModel == upperModelName
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
								 joinKey = nextjoinKey, joinModel = upperModelName))
			} else {
				bMatName <- priorBetween[ which(sameJoinMask) ]
			}

			bMat <- model[[ bMatName ]]
			bMat <- imxGentleResize(bMat, list(rownames(A), upperVars))

			bMat$values[to, fromUpperVar] <- nextvalue
			bMat$free  [to, fromUpperVar] <- nextfree
			bMat$labels[to, fromUpperVar] <- nextlabel
			bMat$ubound[to, fromUpperVar] <- nextubound
			bMat$lbound[to, fromUpperVar] <- nextlbound

			model[[ bMatName ]] <- bMat
			model$expectation$between <- unique(c(priorBetween, bMatName))
			next
		}

		allFromTo <- unique(c(from, to))
		varExist <- allFromTo %in% legalVars
		if(!all(varExist)) {
			fromComponents <- length(strsplit(from, imxSeparatorChar, fixed = TRUE)[[1]])
			if (fromComponents == 2) {
				stop(paste0("In model ", omxQuotes(model@name),
					   ", you tried to add a path from ", omxQuotes(from),
					   ". Did you forget joinKey?"))
			}

			missingVars <- allFromTo[!varExist]
			stop(paste("You need to add",
				   omxQuotes(missingVars),
				   "to either manifestVars or latentVars before you",
				   "can use them in a path."), call. = FALSE)
		}

		if (from == 'one') {
			if (arrows != 1) {
				stop(paste('The means path must be a single-headed arrow\n',
					   'path from "one" to', omxQuotes(to)), call. = FALSE)
			}
		  if (isProductNode[to]) {
		    stop(paste('Cannot change mean of product node', omxQuotes(to),
		               'from the identity value of 1'), call.=FALSE)
		  }
			M@free[1, to] <- nextfree
			M@values[1, to] <- nextvalue
			M@labels[1, to] <- nextlabel
			M@ubound[1, to] <- nextubound
			M@lbound[1, to] <- nextlbound
			next
		}

    if (arrows == 0) {
      # directly modify the unfiltered covariance matrix
      # expectation: data.frame w/ step, rowname, colname ; name of mxMatrix
      # mxMatrix: vector of parameters to plop in
      pair <- c(from,to)[order(match(c(from,to), legalVars))]  # canonical order
      r1 <- data.frame(step=nextstep, from=pair[1], to=pair[2])
      selPlan <- rbind(selPlan, r1)
      oldSelVec <- selVec
      selVec <- mxMatrix('Full', nrow(selPlan), 1)
      if (!is.null(oldSelVec)) selVec[1:nrow(oldSelVec),1] <- oldSelVec
      selVec[nrow(selVec),1]$free <- nextfree
      selVec[nrow(selVec),1]$values <- nextvalue
      selVec[nrow(selVec),1]$labels <- nextlabel
      selVec[nrow(selVec),1]$ubound <- nextubound
      selVec[nrow(selVec),1]$lbound <- nextlbound
      p1 <- order(selPlan$step, selPlan$from, selPlan$to)
      selPlan <- selPlan[p1,]
      selVec <- selVec[p1,]
    } else if (arrows == 1) {
			A_free[to, from] <- nextfree
			A_values[to, from] <- nextvalue
			A_labels[to, from] <- nextlabel
			A_ubound[to, from] <- nextubound
			A_lbound[to, from] <- nextlbound
			# Here, we only check the [to,from] elements (and not the [from,to] elements) of the S_* matrices, 
			# because we assume that 'S' is symmetric, and there's no need to warn about the same thing twice:
			if(S_values[to, from]!=0 || !is.na(S_labels[to, from]) || S_free[to, from]){
				msg <- paste(
					"Looks like there is a pre-existing two-headed path between ",omxQuotes(from)," and ",omxQuotes(to),".\n",
					"That path is now overwritten by a one-headed path from ",omxQuotes(from)," to ",omxQuotes(to),".\n",
					"To retain the two-headed path, either use 'dummy' latent variables, or directly modify the MxModel's 'S' matrix;\n",
					"See the `mxPath()` help page for examples.\n",
					"Be advised, this overwriting behavior may change in the future, so do not write scripts that rely upon it!",
					sep=""
				)
				warning(msg)
			}
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
			if(A_values[to, from]!=0 || !is.na(A_labels[to, from]) || A_free[to, from]){
				msg <- paste(
					"Looks like there is a pre-existing one-headed path from ",omxQuotes(from)," to ",omxQuotes(to),".\n",
					"That path is now overwritten by a two-headed path between ",omxQuotes(from)," and ",omxQuotes(to),".\n",
					"To retain the one-headed path, either use 'dummy' latent variables, or directly modify the MxModel's 'A' matrix;\n",
					"See the `mxPath()` help page for examples.\n",
					"Be advised, this overwriting behavior may change in the future, so do not write scripts that rely upon it!",
					sep=""
				)
				warning(msg)
			}
			A_values[to, from] <- 0
			A_labels[to, from] <- as.character(NA)
			A_free[to, from] <- FALSE
			if(A_values[from, to]!=0 || !is.na(A_labels[from, to]) || A_free[from, to]){
				msg <- paste(
					"Looks like there is a pre-existing one-headed path from ",omxQuotes(to)," to ",omxQuotes(from),".\n",
					"That path is now overwritten by a two-headed path between ",omxQuotes(to)," and ",omxQuotes(from),".\n",
					"To retain the one-headed path, either use 'dummy' latent variables, or directly modify the MxModel's 'A' matrix;\n",
					"See the `mxPath()` help page for examples.\n",
					"Be advised, this overwriting behavior may change in the future, so do not write scripts that rely upon it!",
					sep=""
				)
				warning(msg)
			}
			A_values[from, to] <- 0
			A_labels[from, to] <- as.character(NA)
			A_free[from, to] <- FALSE
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
	model[['A']] <- A
	model[['S']] <- S
	if (!is.null(M)) model[['M']] <- M
  if (!is.null(selVec)) {
    model$expectation$selectionVector <- 'selectionVector'
    model[['selectionVector']] <- selVec
  }
  if (!is.null(selPlan)) model$expectation$selectionPlan <- selPlan
	model
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
	values <- matrix(as.numeric(model$expectation$isProductNode), 1, len)
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
	if (!is.null(model@data)) {
		manifestNames <- observedDataNames(model@data)
		extraVars <- setdiff(model@manifestVars, manifestNames)
    # would need to exclude freq, weight, keys, etc TODO
		## extraData <- setdiff(manifestNames, model@manifestVars)
		## if (length(extraData) > 0) {
		## 	msg <- paste("The observed data contains the variables:",
		## 		omxQuotes(extraData), "that have not been declared in the",
		## 		"manifest variables.")
		## 	stop(msg, call. = FALSE)
		## }
		if (length(extraVars) > 0) {
			msg <- paste("The manifest variables include",
				omxQuotes(extraVars), "that have not been found in the",
				"observed data.")
			stop(msg, call. = FALSE)
		}
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
  # Get the number of manifest variables that are already in the model
  currentManifest <- max(c(0, length(model@manifestVars) - newManifest))
  # Get the number of latent variables that are already in the model
  currentLatent <- max(c(0, length(model@latentVars) - newLatent))
  # The matrix will have the dimensions (n_manifest + n_latent) * (n_manifest + n_latent)
  newSize <- length(model@manifestVars) + length(model@latentVars)
  if (currentManifest > 0) {
    # The structure of the old matrix is (current manifest + current latent) *
    # (current manifest + current latent) 
    # We will keep the curret manifest * manifest part:
    manifestXmanifest <- oldmatrix[1 : currentManifest, 1 : currentManifest, drop = FALSE]
  } else {
    # If there are no manifest variables: Just set up an empty matrix
    manifestXmanifest <- matrix(value, currentManifest, currentManifest)
  }
  if (currentLatent > 0) {
    # If there are latent variables in the model already, we extract the three matrixes
    latentStart <- currentManifest + 1
    latentEnd <- currentManifest + currentLatent
    # (1) manifest * latent
    manifestXlatent <- oldmatrix[1 : currentManifest, latentStart : latentEnd, drop = FALSE]
    # (2) latent * manifest
    latentXmanifest <- oldmatrix[latentStart : latentEnd, 1 : currentManifest, drop = FALSE]
    # (3) latent * latent
    latentXlatent <- oldmatrix[latentStart : latentEnd, latentStart : latentEnd, drop = FALSE]
  } else {
    manifestXlatent <- matrix(value, currentManifest, currentLatent)
    latentXmanifest <- matrix(value, currentLatent, currentManifest)
    latentXlatent <- matrix(value, currentLatent, currentLatent)
  }
  
  newtop <- cbind(
    # Extend the current manifest * manifest matrix by adding a manifest * new manifest
    # matrix
    manifestXmanifest, matrix(value, currentManifest, newManifest),
    # Extend the current manifest * latent by adding a manifest * new latent matrix
    manifestXlatent, matrix(value, currentManifest, newLatent))
  # So far, we have only added columns for the new manifest variables. We also
  # have to add a row for this new manifest variable. This is a 
  # new latent * (all variables) matrix.
  newtop <- rbind(newtop, matrix(value, newManifest, newSize))
  # The lower part of our matrix is dedicated to the latent variables
  newbottom <- cbind(
    # Extend the latent * manifest variables matrix to a latent * (manifest + new manifest)
    # matrix
    latentXmanifest, matrix(value, currentLatent, newManifest),
    # Extend the latent * latent matrix to a latent * (latent + new latent) matrix
    latentXlatent, matrix(value, currentLatent, newLatent))
  # Now, we also add a new latent * (all variables) matrix
  newbottom <- rbind(newbottom, matrix(value, newLatent, newSize))
  
  newmatrix <- rbind(newtop, newbottom)
  return(newmatrix)
}

addVariablesMatrixM <- function(oldmatrix, newLatentValue, newManifestValue, model, newLatent, newManifest) {
  newManifest <- length(newManifest)
  newLatent <- length(newLatent)
  currentManifest <- length(model@manifestVars) - newManifest
  currentLatent <- length(model@latentVars) - newLatent
  if(currentManifest > 0){
    oldManifest <- oldmatrix[1, 1:currentManifest, drop = FALSE]
  }else{
    oldManifest <- matrix(0, 0, 0)
  }
  if(currentLatent > 0){
    oldLatent <- oldmatrix[1, (currentManifest + 1) : (currentLatent + currentManifest), drop = FALSE]
  }else{
    oldLatent <- matrix(0, 0, 0)
  }
  
  values <- c(oldManifest,
              rep.int(newManifestValue, newManifest),
              oldLatent,
              rep.int(newLatentValue, newLatent))
  newmatrix <- matrix(values, 1, length(model@manifestVars) + length(model@latentVars))
  return(newmatrix)
}
