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

observedStatistics <- function(model, data) {
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
		return(length(data@observed))
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
	retval[['AIC']] <- likelihood - 2 * DoF
	retval[['BIC']] <- 0.5 * (likelihood - DoF * log(data@numObs))
	rmseaSquared <- (chi / DoF - 1) / data@numObs
	if (rmseaSquared < 0) {
		retval[['RMSEA']] <- 0
	} else {
		retval[['RMSEA']] <- sqrt(rmseaSquared)
	}
	return(retval)
}

computeOptimizationStatistics <- function(model, matrices, parameters, objective, data) {
	retval <- list()
	if(length(model@output) == 0) { return(retval) }
	ptable <- data.frame()
	estimates <- model@output$estimate
	if (length(estimates) > 0) {
		matrixNames <- names(matrices)
		for(i in 1:length(estimates)) {
			mLocation <- parameters[[i]][[3]][[1]] + 1
			mRow <- parameters[[i]][[3]][[2]] + 1
			mCol <- parameters[[i]][[3]][[3]] + 1
			ptable[i, 'name'] <- names(estimates)[[i]]
			ptable[i, 'matrix'] <- simplifyName(matrixNames[[mLocation]], model@name)
			ptable[i, 'row'] <- mRow
			ptable[i, 'col'] <- mCol
			ptable[i, 'parameter estimate'] <- estimates[[i]]
			ptable[i, 'error estimate'] <- model@output$hessian[i, i]

		}
		retval[['parameters']] <- ptable
	}
	retval[['estimatedParameters']] <- length(estimates)
	retval[['observedStatistics']] <- observedStatistics(model, data)
	retval[['degreesOfFreedom']] <- retval[['observedStatistics']] - retval[['estimatedParameters']]
	retval[['SaturatedLikelihood']] <- model@output$SaturatedLikelihood
	retval[['Minus2LogLikelihood']] <- model@output$Minus2LogLikelihood
	retval <- fitStatistics(model, objective, data, retval)
	return(retval)
}

setMethod("summary", "MxModel",
	function(object, ...) {	
		namespace <- omxGenerateNamespace(object)
		flatModel <- omxFlattenModel(object, namespace)
		matrices <- generateSimpleMatrixList(flatModel)
		parameters <- generateParameterList(flatModel)
		objective <- flatModel@objectives[[omxIdentifier(object@name, 'objective')]]
		data <- flatModel@datasets[[objective@data]]
		retval <- computeOptimizationStatistics(object, matrices, parameters, objective, data)
		if (!is.null(data)) {
			print(summary(data@observed))
			cat('\n')
		}
		if (length(object@output) > 0) {
			message <- npsolMessages[[as.character(object@output$status[[1]])]]
			if (!is.null(message)) {
				cat(message,'\n','\n')
			}
		}
		if (!is.null(retval$parameters)) {
			print(retval$parameters)
			cat('\n')
		}
		cat("Observed statistics: ", retval$observedStatistics, '\n')
		cat("Estimated parameters: ", retval$estimatedParameters, '\n')
		cat("Degrees of freedom: ", retval$degreesOfFreedom, '\n')
		cat("AIC: ", retval$AIC, '\n')
		cat("BIC: ", retval$BIC, '\n')
		cat("adjusted BIC:", '\n')
		cat("RMSEA: ", retval$RMSEA, '\n')
		cat('\n')		
		invisible(retval)
	}
)

