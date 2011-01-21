#
#   Copyright 2007-2010 The OpenMx Project
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

calculateConstraints <- function(model, useSubmodels) {
	constraints <- model@runstate$constraints
	retval <- c()
	if (length(constraints) > 0) {
		retval <- sapply(constraints, calculateConstraintsHelper, model)
	}
	if (useSubmodels && length(model@runstate$independents) > 0) {
		submodelConstraints <- lapply(model@runstate$independents, calculateConstraints, FALSE)
		names(submodelConstraints) <- NULL
		submodelConstraints <- unlist(submodelConstraints)
		retval <- c(retval, submodelConstraints)
	}
	return(retval)
}

calculateConstraintsHelper <- function(constraint, model) {
	if (constraint@relation == "==") {
		leftHandSide <- constraint@formula[[2]]
		value <- eval(substitute(mxEval(x, model, compute=TRUE),
			list(x = leftHandSide)))
		value <- as.matrix(value)
		return(nrow(value) * ncol(value))
	} else {
		return(0)
	}
}

observedStatisticsHelper <- function(objective, datalist, historySet) {
	if ('numStats' %in% slotNames(objective)) {
		if (!is.na(objective@numStats)) {
			return(list(objective@numStats, historySet))
		}
	}
	if (is.na(objective@data)) {
		return(list(0, historySet))
	}
	data <- datalist[[objective@data + 1]]
	if (data@type == 'cov' || data@type == 'sscp') {
		if (data@name %in% historySet) {
			return (list(0, historySet))
		}
		n <- nrow(data@observed)
		dof <- n * (n + 1) / 2
		if (!single.na(data@means)) {
			dof <- dof + length(data@means)
		}
		historySet <- append(data, historySet)		
	} else if (data@type == 'cor') {
		if (data@name %in% historySet) {
			return (list(0, historySet))
		}
		n <- nrow(data@observed)
		dof <- n * (n - 1) / 2
		if (!single.na(data@means)) {
			dof <- dof + length(data@means) 
		}
		historySet <- append(data, historySet)
	} else {
		dof <- 0
		observed <- data@observed
		for (i in 1:ncol(observed)) {
			colname <- colnames(observed)[[i]]
			fullname <- paste(data@name, colname, sep='.')
			if ((colname %in% objective@dims) && !(fullname %in% historySet)) {
				dof <- dof + sum(!is.na(observed[,i]))
				historySet <- append(fullname, historySet)
			}
		}
	}
	return(list(dof, historySet))
}

observedStatistics <- function(model, useSubmodels, constraintOffset) {
	datalist <- model@runstate$datalist
	objectives <- model@runstate$objectives
	retval <- constraintOffset
	if (length(objectives) > 0) {
		historySet <- character()
		for(i in 1:length(objectives)) {
			result <- observedStatisticsHelper(objectives[[i]], datalist, historySet)
			retval <- retval + result[[1]]
			historySet <- result[[2]]
		}
	}
	if (useSubmodels && length(model@runstate$independents) > 0) {
		submodelStatistics <- sapply(model@runstate$independents, observedStatistics, FALSE, 0)
		retval <- retval + sum(submodelStatistics)
	}
	return(retval)
}

objectiveNumberObservations <- function(objective) {
	if (("numObs" %in% slotNames(objective)) && !single.na(objective@numObs)) {
		return(objective@numObs)
	} else {
		return(0)
	}
}

numberObservations <- function(datalist, objectives) {
	dataObservations <- sapply(datalist, slot, name = "numObs")
	objectiveObservations <- sapply(objectives, objectiveNumberObservations)
	return(sum(as.numeric(dataObservations), as.numeric(objectiveObservations)))
}

computeFValue <- function(datalist, likelihood, chi) {
	if(length(datalist) == 0) return(NA)
	if(all(sapply(datalist, function(x) 
		{x@type == 'raw'}))) return(likelihood)
	if(all(sapply(datalist, function(x) 
		{x@type == 'cov'}))) return(chi)
	return(NA)
}

fitStatistics <- function(model, useSubmodels, retval) {
	datalist <- model@runstate$datalist
	likelihood <- retval[['Minus2LogLikelihood']]
	saturated <- retval[['SaturatedLikelihood']]
	chi <- likelihood - saturated
	DoF <- retval$degreesOfFreedom
	Fvalue <- computeFValue(datalist, likelihood, chi)
	retval[['Chi']] <- chi
	retval[['p']] <- suppressWarnings(pchisq(chi, DoF, lower.tail = FALSE))
	retval[['AIC.Mx']] <- Fvalue - 2 * DoF
	retval[['BIC.Mx']] <- 0.5 * (Fvalue - DoF * log(retval[['numObs']]))
	rmseaSquared <- (chi / DoF - 1) / retval[['numObs']]
	if (length(rmseaSquared) == 0 || is.na(rmseaSquared) || 
		is.nan(rmseaSquared) || (rmseaSquared < 0)) {
		retval[['RMSEA']] <- NA
	} else {
		retval[['RMSEA']] <- sqrt(rmseaSquared)
	}
	return(retval)
}

parameterList <- function(model, useSubmodels) {
	if (useSubmodels && length(model@runstate$independents) > 0) {
		ptable <- parameterListHelper(model, TRUE)
		submodelParameters <- lapply(model@runstate$independents, parameterListHelper, TRUE)
		ptable <- Reduce(rbind, submodelParameters, ptable)
	} else {
		ptable <- parameterListHelper(model, FALSE)
	}
	return(ptable)
}

parameterListHelper <- function(model, modelName) {
	ptable <- data.frame()
	if(length(model@output) == 0) { return(ptable) }
	estimates <- model@output$estimate
    if (!is.null(model@output$standardErrors) && 
 		length(model@output$standardErrors) == length(estimates)) { 
 		errorEstimates <- model@output$standardErrors 
	} else { 
		errorEstimates <- rep.int(as.numeric(NA), length(estimates)) 
	}
	matrices <- model@runstate$matrices
	parameters <- model@runstate$parameters
	if (length(estimates) > 0) {
		matrixNames <- names(matrices)
		for(i in 1:length(estimates)) {
			mLocation <- parameters[[i]][[3]][[1]] + 1
			mRow <- parameters[[i]][[3]][[2]] + 1
			mCol <- parameters[[i]][[3]][[3]] + 1
			lbound <- parameters[[i]][[1]]
			ubound <- parameters[[i]][[2]]
			aMatrix <- matrices[[mLocation]][[1]]
			if (getOption('mxShowDimnames')) {
				if (!is.null(rownames(aMatrix))) {
					mRow <- rownames(aMatrix)[[mRow]]
				}
				if (!is.null(colnames(aMatrix))) {
					mCol <- colnames(aMatrix)[[mCol]]
				}
			}
			if (modelName) {
				ptable[i, 'model'] <- model@name
			}
			ptable[i, 'name'] <- names(estimates)[[i]]
			ptable[i, 'matrix'] <- simplifyName(matrixNames[[mLocation]], model@name)
			ptable[i, 'row'] <- mRow
			ptable[i, 'col'] <- mCol
			ptable[i, 'Estimate'] <- estimates[[i]]
			ptable[i, 'Std.Error'] <- errorEstimates[[i]]
			ptable[i, 'lbound'] <- lbound
			ptable[i, 'ubound'] <- ubound
		}
	}
	return(ptable)
}

computeOptimizationStatistics <- function(model, numStats, useSubmodels, retval) {
	estimates <- model@output$estimate
	retval[['constraints']] <- calculateConstraints(model, useSubmodels)
	retval[['estimatedParameters']] <- nrow(retval$parameters)
	if (is.null(numStats)) {
		retval[['observedStatistics']] <- observedStatistics(model, useSubmodels, sum(retval$constraints))
	} else {
		retval[['observedStatistics']] <- numStats
	}
	retval[['degreesOfFreedom']] <- retval$observedStatistics - retval$estimatedParameters
	retval <- fitStatistics(model, useSubmodels, retval)
	return(retval)
}

print.summary.mxmodel <- function(x,...) {
	if (length(x$dataSummary) > 0) {
		cat("data:\n")
		print(x$dataSummary)
	}
	if (!is.null(x$npsolMessage)) {
		cat(x$npsolMessage,'\n','\n')
	}
	if (length(x$parameters) > 0) {
		cat("free parameters:\n")
		params <- x$parameters
		params$lbound[is.na(params$lbound)] <- ""
		params$ubound[is.na(params$ubound)] <- ""
		print(params)
		cat('\n')
	}
	if (!is.null(x$CI) && length(x$CI) > 0) {
		cat("confidence intervals:\n")
		print(x$CI)
		cat('\n')
	}
	cat("observed statistics: ", x$observedStatistics, '\n')
	constraints <- x$constraints
	if(length(constraints) > 0) {
		for(i in 1:length(constraints)) {
			name <- names(constraints)[[i]]
			if (constraints[[i]] == 1) plural <- ''
			else plural <- 's'
			cat("Constraint", omxQuotes(simplifyName(name, x$modelName)), "contributes",
				constraints[[i]], paste("observed statistic", plural, '.', sep=''), "\n")
		}
	}
	cat("estimated parameters: ", x$estimatedParameters, '\n')
	cat("degrees of freedom: ", x$degreesOfFreedom, '\n')
	cat("-2 log likelihood: ", x$Minus2LogLikelihood, '\n')
	cat("saturated -2 log likelihood: ", x$SaturatedLikelihood, '\n')
	cat("number of observations: ", x$numObs, '\n')
	cat("chi-square: ", x$Chi, '\n')
	cat("p: ", x$p, '\n')
	cat("AIC (Mx): ", x$AIC.Mx, '\n')
	cat("BIC (Mx): ", x$BIC.Mx, '\n')
	cat("adjusted BIC:", '\n')
	cat("RMSEA: ", x$RMSEA, '\n')
	cat("timestamp:", format(x$timestamp), '\n')
	cat("frontend time:", format(x$frontendTime), '\n')
	cat("backend time:", format(x$backendTime), '\n')
	cat("independent submodels time:", format(x$independentTime), '\n')
	cat("wall clock time:", format(x$wallTime), '\n')
	cat("cpu time:", format(x$cpuTime), '\n')
	cat("openmx version number:", format(x$mxVersion), '\n')
	cat('\n')
}

setLikelihoods <- function(model, saturatedLikelihood, retval) {
	if(is.null(saturatedLikelihood)) {
		retval$SaturatedLikelihood <- model@output$SaturatedLikelihood
	} else {
		retval$SaturatedLikelihood <- saturatedLikelihood		
	}
	retval$Minus2LogLikelihood <- model@output$Minus2LogLikelihood
	if (is.null(retval$SaturatedLikelihood)) {
		retval$SaturatedLikelihood <- NA
	}
	if (is.null(retval$Minus2LogLikelihood)) {
		retval$Minus2LogLikelihood <- NA
	}
	return(retval)
}

setNumberObservations <- function(numObs, datalist, objectives, retval) {
	if(is.null(numObs)) {
		retval$numObs <- numberObservations(datalist, objectives)
	} else {
		retval$numObs <- numObs
	}
	return(retval)
}

summarizeDataObject <- function(dataObject) {
	if (dataObject@type != "raw") {
		result <- list()
		result[[dataObject@type]] <- dataObject@observed
		if (!single.na(dataObject@means)) {
			result[['means']] <- dataObject@means
		}
		return(result)
	} else {
		return(summary(dataObject@observed))
	}
}

generateDataSummary <- function(model, useSubmodels) {
	datalist <- model@runstate$datalist
	retval <- lapply(datalist, summarizeDataObject)
	if (useSubmodels && length(model@runstate$independents) > 0) {
		submodelSummary <- lapply(model@runstate$independents, generateDataSummary, FALSE)
		names(submodelSummary) <- NULL
		submodelSummary <- unlist(submodelSummary, recursive = FALSE)
		retval <- c(retval, submodelSummary)
	}
	return(retval)
}

imxEvalByName <- function(name, model, compute=FALSE, show=FALSE) {
   if ((length(name) != 1) || typeof(name) != "character") {
      stop("'name' argument must be a character argument")
   }
   if (!is(model, "MxModel")) {
      stop("'model' argument must be a MxModel object")
   }
   if (hasSquareBrackets(name)) {
      components <- splitSubstitution(name)
      eval(substitute(mxEval(x[y,z], model, compute, show),
         list(x = as.name(components[[1]]), 
            y = parse(text=components[[2]])[[1]],
            z = parse(text=components[[3]])[[1]])))
   } else {
      eval(substitute(mxEval(x, model, compute, show),
         list(x = as.name(name))))
   }
}

generateConfidenceIntervalTable <- function(model) {
	base <- model@output$confidenceIntervals
	if (length(base) == 0) return(matrix(0, 0, 3))
	entities <- rownames(base)
	estimates <- sapply(entities, imxEvalByName, model, compute=TRUE, show=FALSE)
	retval <- cbind(base[, 'lbound'], estimates, base[, 'ubound'])
	rownames(retval) <- entities
	colnames(retval) <- c('lbound', 'estimate', 'ubound')
	return(retval)
}

translateSaturatedLikelihood <- function(input) {
	if (is.null(input)) {
		return(input)
	} else if (is.numeric(input)) {
		return(input)
	} else if (is(input, "MxModel")) {
		if (is.null(input@objective)) {
			stop(paste("Saturated model passed",
				"to summary function does not",
				"have top-level objective function in",
				deparse(width.cutoff = 400L, sys.call(-1))), call. = FALSE)
		}
		if (length(input@objective@result) != 1) {
			stop(paste("Saturated model passed to summary",
				"function does not have a 1x1 matrix",
				"result in objective function in",
				deparse(width.cutoff = 400L, sys.call(-1))), call. = FALSE)
		}
		return(input@objective@result[1,1])
	} else {
		stop(paste("Illegal argument passed to",
			"'SaturatedLikelihood' argument",
			"of summary function in",
			deparse(width.cutoff = 400L, sys.call(-1))), call. = FALSE)
	}
}

setMethod("summary", "MxModel",
	function(object, ...) {
		model <- object
		dotArguments <- list(...)
		saturatedLikelihood <- translateSaturatedLikelihood(dotArguments$SaturatedLikelihood)
		numObs <- dotArguments$numObs
		numStats <- dotArguments$numStats
		useSubmodels <- dotArguments$indep
		if (is.null(useSubmodels)) { useSubmodels <- TRUE }
		retval <- list()
		retval$parameters <- parameterList(model, useSubmodels)
		retval <- setLikelihoods(model, saturatedLikelihood, retval)
		retval <- setNumberObservations(numObs, model@runstate$datalist, model@runstate$objectives, retval)
		retval <- computeOptimizationStatistics(model, numStats, useSubmodels, retval)
		retval$dataSummary <- generateDataSummary(model, useSubmodels)
		retval$CI <- generateConfidenceIntervalTable(model)
		retval$CIcodes <- model@output$confidenceIntervalCodes
		if (!is.null(model@output$status)) {
			message <- npsolMessages[[as.character(model@output$status[[1]])]]
			retval[['npsolMessage']] <- message
		}
		retval$timestamp <- model@output$timestamp
		retval$frontendTime <- model@output$frontendTime
		retval$backendTime <- model@output$backendTime
		retval$independentTime <- model@output$independentTime
		retval$wallTime <- model@output$wallTime
		retval$cpuTime <- model@output$cpuTime
		retval$mxVersion <- model@output$mxVersion
		retval$modelName <- model@name
		class(retval) <- "summary.mxmodel"
		return(retval)
	}
)
