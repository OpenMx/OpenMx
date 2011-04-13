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

modelIsHollow <- function(model) {
	return(is.null(model@objective) && 
		length(model@matrices) == 0 && 
		length(model@algebras) == 0 &&
		length(omxDependentModels(model)) == 0)
}

processHollowModel <- function(model, independents, dataList, frontendStart, indepElapsed) {
	independents <- lapply(independents, undoDataShare, dataList)
	model <- omxReplaceModels(model, independents)
	model <- undoDataShare(model, dataList)
	frontendStop <- Sys.time()
	frontendElapsed <- (frontendStop - frontendStart) - indepElapsed
	model@output <- calculateTiming(model@output, frontendElapsed, 
		0, indepElapsed, frontendStop, independents)
	model@output$mxVersion <- mxVersion()
	model@runstate$independents <- independents
	return(model)
}


processOptimizerOutput <- function(suppressWarnings, flatModel, matrixNames, 
		algebraNames, parameterNames, intervalNames, unsafe, output) {
	output$mxVersion <- mxVersion()
	if (length(output$estimate) == length(parameterNames)) {
		names(output$estimate) <- parameterNames
	}
	if (length(output$gradient) == length(parameterNames)) {
		names(output$gradient) <- parameterNames
	}
	output$estimatedHessian <- t(output$hessianCholesky) %*% output$hessianCholesky
	if (nrow(output$estimatedHessian) == length(parameterNames) &&
		ncol(output$estimatedHessian) == length(parameterNames)) {
		dimnames(output$estimatedHessian) <- list(parameterNames, parameterNames)
	}
	if (!is.null(output$calculatedHessian) && 
		nrow(output$calculatedHessian) == length(parameterNames) &&
		ncol(output$calculatedHessian) == length(parameterNames)) {
		dimnames(output$calculatedHessian) <- list(parameterNames, parameterNames)
	}
	if (length(output$matrices) == length(matrixNames)) {
		names(output$matrices) <- matrixNames
	}
	if (length(output$algebras) == length(algebraNames)) {
		names(output$algebras) <- algebraNames
	}
	if (length(output$confidenceIntervals) > 0) {
		dimnames(output$confidenceIntervals) <- list(intervalNames, c('lbound', 'ubound'))
		dimnames(output$confidenceIntervalCodes) <- list(intervalNames, c('lbound', 'ubound'))
	}
	if (output$status[[1]] > 0 && !suppressWarnings) {
		npsolWarnings(paste("In model", omxQuotes(flatModel@name)), output$status[[1]])
	}
	if (output$status[[1]] < 0 || output$status[[2]] < 0) {
		if (unsafe) {
			warning(paste("The job for model", 
			omxQuotes(flatModel@name),
			"exited abnormally with the error message:",
			output$status[[3]]), call. = FALSE)
		} else {
			stop(paste("The job for model",
			omxQuotes(flatModel@name),
			"exited abnormally with the error message:",
			output$status[[3]]), call. = FALSE)
		}
	}
	return(output)
}

populateRunStateInformation <- function(model, parameters, matrices, 
		objectives, datalist, constraints, independents, defvars) {
	model@runstate$parameters <- parameters
	model@runstate$matrices <- matrices
	model@runstate$objectives <- objectives
	model@runstate$datalist <- datalist
	model@runstate$constraints <- constraints
	model@runstate$independents <- independents
	model@runstate$defvars <- names(defvars)
	return(model)
}

calculateTiming <- function(output, frontend,
	backend, indep, timestamp, independents) {
	output$frontendTime <- frontend
	output$backendTime <- backend
	output$independentTime <- indep
	output$wallTime <- frontend + backend + indep
	output$timestamp <- timestamp
	if("package:snowfall" %in% search() && length(independents) > 0) {
		output$cpuTime <- frontend + backend + sum(sapply(independents,
			function(x) { x@output$cpuTime }))
	} else {
		output$cpuTime <- output$wallTime
	}
	return(output)
}

npsolMessages <- list('1' = paste('The final iterate satisfies',
		'the optimality conditions to the accuracy requested,',
		'but the sequence of iterates has not yet converged.',
		'NPSOL was terminated because no further improvement',
		'could be made in the merit function (Mx status GREEN).'),
		'2' = paste('The linear constraints and bounds could not be satisfied.',
		'The problem has no feasible solution.'),
		'3' = paste('The nonlinear constraints and bounds could not be satisfied.',
		'The problem may have no feasible solution.'),
		'4' = 'The major iteration limit was reached (Mx status BLUE).',
		'6' = paste('The model does not satisfy the first-order optimality conditions',
		'to the required accuracy, and no improved point for the',
		'merit function could be found during the final linesearch (Mx status RED)'),
		'7' = paste('The function derivates returned by funcon or funobj',
		'appear to be incorrect.'),
		'9' = 'An input parameter was invalid')

npsolWarnings <- function(prefix, status) {
	message <- npsolMessages[[as.character(status)]]
	if(!is.null(message)) {
		warning(paste(prefix, 
			"NPSOL returned a non-zero status code", 
			paste(status, '.', sep = ''), message), call. = FALSE)
	}
}
