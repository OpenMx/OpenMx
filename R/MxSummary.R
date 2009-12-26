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

observedStatistics <- function(model, data, flatModel) {
	if (is.null(data)) {
		return(0)
	}
	if (data@type == 'cov' || data@type == 'sscp') {
		n <- nrow(data@observed)
		dof <- n * (n + 1) / 2 + length(flatModel@constraints)
		if (!single.na(data@means)) {
			dof <- dof + length(data@means)
		}
		return(dof)
	} else if (data@type == 'cor') {
		n <- nrow(data@observed)
		dof <- n * (n - 1) / 2 + length(flatModel@constraints)
		if (!single.na(data@means)) {
			dof <- dof + length(data@means) 
		}
		return(dof)
	} else {
		return(sum(!is.na(data@observed)) + length(flatModel@constraints))
	}
}

fitStatistics <- function(model, objective, data, retval) {
	likelihood <- model@output$Minus2LogLikelihood
	saturated <- model@output$SaturatedLikelihood
	chi <- likelihood - saturated
	DoF <- retval$degreesOfFreedom
	if (is.null(objective) || is(objective, "MxAlgebraObjective")) {
		return(retval)
	} else if (is.null(likelihood)) {
		return(retval)
	}
	if (is.null(data)) {
		Fvalue <- NA
	} else if (data@type == 'raw') {
		Fvalue <- likelihood
	} else if (data@type == 'cov') {
		Fvalue <- chi
	} else {
		Fvalue <- NA
	}
	if (!is.null(data)) {
		retval[['Chi']] <- chi
		retval[['p']] <- pchisq(chi, DoF, lower.tail = FALSE)
		retval[['AIC.Mx']] <- Fvalue - 2 * DoF
		retval[['BIC.Mx']] <- 0.5 * (Fvalue - DoF * log(data@numObs))
		rmseaSquared <- (chi / DoF - 1) / data@numObs
		if (length(rmseaSquared) == 0 || is.nan(rmseaSquared) || (rmseaSquared < 0)) {
			retval[['RMSEA']] <- 0
		} else {
			retval[['RMSEA']] <- sqrt(rmseaSquared)
		}
	}
	return(retval)
}

computeOptimizationStatistics <- function(model, matrices, parameters, objective, data, flatModel) {
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
	retval[['estimatedParameters']] <- length(estimates)
	retval[['observedStatistics']] <- observedStatistics(model, data, flatModel)
	retval[['degreesOfFreedom']] <- retval[['observedStatistics']] - retval[['estimatedParameters']]
	retval[['SaturatedLikelihood']] <- model@output$SaturatedLikelihood
	retval[['Minus2LogLikelihood']] <- model@output$Minus2LogLikelihood
	retval <- fitStatistics(model, objective, data, retval)
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
	cat("Chi-Square: ", x$Chi, '\n')
	cat("p: ", x$p, '\n')
	cat("AIC (Mx): ", x$AIC.Mx, '\n')
	cat("BIC (Mx): ", x$BIC.Mx, '\n')
	cat("adjusted BIC:", '\n')
	cat("RMSEA: ", x$RMSEA, '\n')
	cat('\n')
}

setMethod("summary", "MxModel",
	function(object, ...) {	
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
		retval <- computeOptimizationStatistics(object, matrices, parameters, objective, data, flatModel)
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