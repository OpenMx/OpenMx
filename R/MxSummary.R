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

dataMatchesHistory <- function(data, historySet) {
	return(any(sapply(historySet, identical, data)))
}

observedStatisticsHelper <- function(data, historySet) {
	if (is.null(data)) {
		return(list(0, historySet))
	}
	if (data@type == 'cov' || data@type == 'sscp') {
		if (dataMatchesHistory(data, historySet)) {
			return (list(0, historySet))
		}
		n <- nrow(data@observed)
		dof <- n * (n + 1) / 2
		if (!single.na(data@means)) {
			dof <- dof + length(data@means)
		}
		historySet <- append(data, historySet)		
	} else if (data@type == 'cor') {
		if (dataMatchesHistory(data, historySet)) {
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
		for (i in 1:ncol(data@observed)) {
			if (!dataMatchesHistory(data@observed[,i], historySet)) {
				dof <- dof + sum(!is.na(data@observed[,i]))
				historySet <- append(data@observed[,i], historySet)
			}
		}
	}
	return(list(dof, historySet))
}

observedStatisticsSingleModel <- function(data) {
	if (is.null(data)) {
		return(0)
	}
	if (data@type == 'cov' || data@type == 'sscp') {
		n <- nrow(data@observed)
		dof <- n * (n + 1) / 2
		if (!single.na(data@means)) {
			dof <- dof + length(data@means)
		}
		return(dof)
	} else if (data@type == 'cor') {
		n <- nrow(data@observed)
		dof <- n * (n - 1) / 2
		if (!single.na(data@means)) {
			dof <- dof + length(data@means) 
		}
		return(dof)
	} else {
		return(sum(!is.na(data@observed)))
	}
}

observedStatistics <- function(flatModel) {
	datasets <- flatModel@datasets
	if (length(datasets) == 0) {
		return(length(flatModel@constraints))
	} else if (length(datasets) == 1) {
		return(observedStatisticsSingleModel(datasets[[1]]) + 
			length(flatModel@constraints))
	}
	historySet <- list()
	retval <- 0
	for(i in 1:length(datasets)) {
		result <- observedStatisticsHelper(datasets[[i]], historySet)
		retval <- retval + result[[1]]
		historySet <- result[[2]]
	}
	retval <- retval + length(flatModel@constraints)
	return(retval)
}

numberObservations <- function(flatModel) {
	obs <- sapply(flatModel@datasets, 
		function(x) { x@numObs })
	return(sum(as.numeric(obs)))
}

computeFValue <- function(flatModel, likelihood, chi) {
	if(length(flatModel@datasets) == 0) return(NA)
	if(all(sapply(flatModel@datasets, function(x) 
		{x@type == 'raw'}))) return(likelihood)
	if(all(sapply(flatModel@datasets, function(x) 
		{x@type == 'cov'}))) return(chi)
	return(NA)
}

fitStatistics <- function(flatModel, retval) {
	likelihood <- retval[['Minus2LogLikelihood']]
	saturated <- retval[['SaturatedLikelihood']]
	chi <- likelihood - saturated
	DoF <- retval$degreesOfFreedom
	Fvalue <- computeFValue(flatModel, likelihood, chi)
	retval[['Chi']] <- chi
	retval[['p']] <- pchisq(chi, DoF, lower.tail = FALSE)
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

parameterList <- function(model, matrices, parameters) {
	retval <- list()
	if(length(model@output) == 0) { return(retval) }
	ptable <- data.frame()
	estimates <- model@output$estimate
	if (length(model@output$hessian) > 0) {
		errorEstimates <- sqrt(diag(solve(model@output$hessian)))
	} else {
		errorEstimates <- rep.int(NA, length(estimates))
	}
	if (length(estimates) > 0) {
		matrixNames <- names(matrices)
		for(i in 1:length(estimates)) {
			mLocation <- parameters[[i]][[3]][[1]] + 1
			mRow <- parameters[[i]][[3]][[2]] + 1
			mCol <- parameters[[i]][[3]][[3]] + 1
			aMatrix <- matrices[[mLocation]][[1]]
			if (getOption('mxShowDimnames')) {
				if (!is.null(rownames(aMatrix))) {
					mRow <- rownames(aMatrix)[[mRow]]
				}
				if (!is.null(colnames(aMatrix))) {
					mCol <- colnames(aMatrix)[[mCol]]
				}
			}
			ptable[i, 'name'] <- names(estimates)[[i]]
			ptable[i, 'matrix'] <- simplifyName(matrixNames[[mLocation]], model@name)
			ptable[i, 'row'] <- mRow
			ptable[i, 'col'] <- mCol
			ptable[i, 'Estimate'] <- estimates[[i]]
			ptable[i, 'Std.Error'] <- errorEstimates[[i]]

		}
		retval[['parameters']] <- ptable
	}
	return(retval)
}

computeOptimizationStatistics <- function(model, flatModel, retval) {
	estimates <- model@output$estimate
	retval[['estimatedParameters']] <- length(estimates)
	retval[['observedStatistics']] <- observedStatistics(flatModel)
	retval[['degreesOfFreedom']] <- retval[['observedStatistics']] - retval[['estimatedParameters']]
	retval <- fitStatistics(flatModel, retval)
	return(retval)
}

print.summary.mxmodel <- function(x,...) {
	if (!is.null(x$dataSummary)) {
		print(x$dataSummary)
		cat('\n')
	}
	if (!is.null(x$npsolMessage)) {
		cat(x$npsolMessage,'\n','\n')
	}
	if (!is.null(x$parameters)) {
		print(x$parameters)
		cat('\n')
	}
	cat("Observed statistics: ", x$observedStatistics, '\n')
	cat("Estimated parameters: ", x$estimatedParameters, '\n')
	cat("Degrees of freedom: ", x$degreesOfFreedom, '\n')
	cat("-2 log likelihood: ", x$Minus2LogLikelihood, '\n')
	cat("Saturated -2 log likelihood: ", x$SaturatedLikelihood, '\n')
	cat("numObs: ", x$numObs, '\n')
	cat("Chi-Square: ", x$Chi, '\n')
	cat("p: ", x$p, '\n')
	cat("AIC (Mx): ", x$AIC.Mx, '\n')
	cat("BIC (Mx): ", x$BIC.Mx, '\n')
	cat("adjusted BIC:", '\n')
	cat("RMSEA: ", x$RMSEA, '\n')
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

setNumberObservations <- function(numObs, flatModel, retval) {
	if(is.null(numObs)) {
		retval$numObs <- numberObservations(flatModel)
	} else {
		retval$numObs <- numObs
	}
	return(retval)
}

setMethod("summary", "MxModel",
	function(object, ...) {
		saturatedLikelihood <- match.call()$SaturatedLikelihood
		numObs <- match.call()$numObs
		object <- convertSquareBracketLabels(object)
		namespace <- omxGenerateNamespace(object)
		flatModel <- omxFlattenModel(object, namespace)
		matrices <- generateMatrixList(flatModel)
		parameters <- generateParameterList(flatModel)
		objective <- flatModel@objectives[[omxIdentifier(object@name, 'objective')]]
		if (!is.null(objective)) {
			data <- flatModel@datasets[[objective@data]]
		} else {
			data <- NULL
		}
		retval <- parameterList(object, matrices, parameters)
		retval <- setLikelihoods(object, saturatedLikelihood, retval)
		retval <- setNumberObservations(numObs, flatModel, retval)
		retval <- computeOptimizationStatistics(object, flatModel, retval)
		if (!is.null(data)) {
			retval[['dataSummary']] <- summary(data@observed)
		}
		if (length(object@output) > 0) {
			message <- npsolMessages[[as.character(object@output$status[[1]])]]
			retval[['npsolMessage']] <- message
		}
		class(retval) <- "summary.mxmodel"
		return(retval)
	}
)