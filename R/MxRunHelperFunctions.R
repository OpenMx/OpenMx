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

modelIsHollow <- function(model) {
	if (!isHollow(model)) {
		return(FALSE)
	}
	submodels <- imxDependentModels(model)
	if (length(submodels) == 0) return(TRUE)
	children <- sapply(submodels, modelIsHollow)
	return(all(children))	
}

isHollow <- function(model) {
	return(is.null(model@fitfunction) && 
		length(model@matrices) == 0 && 
		length(model@algebras) == 0)
}

checkResetData <- function(model) {
	if (model@.resetdata) {
		model@data <- NULL
		model@.resetdata <- FALSE
	}
	return(model)
}

processParentData <- function(model, parentData) {
	if (is.null(model@data) && !is.null(parentData)) {
		model@data <- parentData
		model@.resetdata <- TRUE
	}
	return(model)
}

processHollowModel <- function(model, independents, frontendStart, indepElapsed) {
    independents <- lapply(independents, checkResetData)
	model <- imxReplaceModels(model, independents)
	frontendStop <- Sys.time()
	frontendElapsed <- (frontendStop - frontendStart) - indepElapsed
	model@output <- calculateTiming(model@output, frontendElapsed, 
		0, indepElapsed, frontendStop, independents)
	model@output$mxVersion <- mxVersion(verbose=FALSE)
	model@runstate$independents <- independents
	return(model)
}


nameOptimizerOutput <- function(suppressWarnings, flatModel, matrixNames, 
		algebraNames, parameterNames, output) {
	output$mxVersion <- mxVersion(verbose=FALSE)
	if (length(output$estimate) == length(parameterNames)) {
		names(output$estimate) <- parameterNames
	}
	if (length(output$gradient) == length(parameterNames)) {
		names(output$gradient) <- parameterNames
	}
	for (deriv in c("calculatedHessian", "hessian", "ihessian")) {
		mat <- output[[deriv]]
		if (!is.null(mat) &&
		    nrow(mat) == length(parameterNames) &&
		    ncol(mat) == length(parameterNames)) {
			dimnames(mat) <- list(parameterNames, parameterNames)
			output[[deriv]] <- mat
		}
	}
	if (length(output$matrices) == length(matrixNames)) {
		names(output$matrices) <- matrixNames
	}
	if (length(output$algebras) == length(algebraNames)) {
		names(output$algebras) <- algebraNames
	}
	if (length(output$standardErrors) == length(parameterNames)) {
		rownames(output$standardErrors) <- parameterNames
	}
	return(output)
}

processErrorConditions <- function(model, unsafe, suppressWarnings) {
	output <- model@output
	if (!is.null(output$status$code) && !suppressWarnings) {
		npsolWarnings(paste("In model", omxQuotes(model@name)), output$status$code)
	}
    if (!is.null(output$error)) {
    	if (unsafe) {
    		warning(paste("The job for model", omxQuotes(model@name),
	            "exited abnormally with the error message:",
	            output$error), call. = FALSE)
    	} else {
	        stop(paste("The job for model", omxQuotes(model@name),
	            "exited abnormally with the error message:",
	            output$error), call. = FALSE)
	    }
    }
}

getCpuTime <- function(model) {
	return(model@output$cpuTime)
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
			getCpuTime))
	} else {
		output$cpuTime <- output$wallTime
	}
	return(output)
}

optimizerMessages <- list('1' = paste('The final iterate satisfies',
		'the optimality conditions to the accuracy requested,',
		'but the sequence of iterates has not yet converged.',
		'Optimizer was terminated because no further improvement',
		'could be made in the merit function (Mx status GREEN).'),
		'2' = paste('The linear constraints and bounds could not be satisfied.',
		'The problem has no feasible solution.'),
		'3' = paste('The nonlinear constraints and bounds could not be satisfied.',
		'The problem may have no feasible solution.'),
		'4' = 'The major iteration limit was reached (Mx status BLUE).',
		'5' = 'The Hessian at the solution does not appear to be convex (Mx status RED).',
		'6' = paste('The model does not satisfy the first-order optimality conditions',
		'to the required accuracy, and no improved point for the',
		'merit function could be found during the final linesearch (Mx status RED)'),
		'7' = paste('The function derivates returned by funcon or funobj',
		'appear to be incorrect.'),
		'9' = 'An input parameter was invalid',
		      '10' = 'Starting values are not feasible. Consider mxTryHard()')

npsolWarnings <- function(prefix, status) {
	message <- optimizerMessages[[as.character(status)]]
	if(!is.null(message)) {
		warning(paste(prefix, 
			"Optimizer returned a non-zero status code",
			paste(status, '.', sep = ''), message), call. = FALSE)
	}
}
