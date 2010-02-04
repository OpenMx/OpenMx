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

mxRun <- function(model, silent = FALSE) {
	frontendStart <- Sys.time()
	if(!silent) cat("Running", model@name, "\n")
	namespace <- omxGenerateNamespace(model)
	omxCheckNamespace(model, namespace)
	omxCheckMatrices(model)
	omxVerifyModel(model)
	model <- translateObjectives(model, namespace)
	# Regenerate the namespace
	namespace <- omxGenerateNamespace(model)
	dshare <- shareData(model)
	independents <- getAllIndependents(dshare)
	independents <- omxLapply(independents, mxRun, silent)
	frozen <- omxLapply(independents, omxFreezeModel)
	model <- omxReplaceModels(model, frozen)
	if(is.null(model@objective) && 
			length(model@matrices) == 0 && 
			length(model@algebras) == 0 &&
			length(omxDependentModels(model)) == 0) {
		frontendStop1 <- Sys.time()
		backendStop <- frontendStop1
		frontendStop2 <- frontendStop1
		model@output <- calculateTiming(model@output, frontendStart,
			frontendStop1, backendStop, frontendStop2)		
		return(model)
	}
	model <- convertSquareBracketLabels(model)
	flatModel <- omxFlattenModel(model, namespace)
	data <- convertDatasets(flatModel, model)
	freeFixedValues <- omxCheckVariables(flatModel, namespace)
	flatModel <- populateDefInitialValues(flatModel)
	oldFlatModel <- flatModel
	flatModel <- convertAlgebras(flatModel, list(startvals=freeFixedValues, 
		values=namespace$values, parameters=namespace$parameters))
	cycleDetection(flatModel)
	flatModel <- checkEvaluation(model, flatModel, oldFlatModel)
	parameters <- generateParameterList(flatModel)
	matrices <- generateMatrixList(flatModel)
	algebras <- generateAlgebraList(flatModel)
	startVals <- generateValueList(matrices, parameters)
	objectives <- convertObjectives(flatModel, model)
	algebras <- append(algebras, objectives)
	constraints <- convertConstraints(flatModel)
	state <- c()
	objective <- getObjectiveIndex(flatModel)
	options <- generateOptionsList(model@options)
	frontendStop1 <- Sys.time()
	output <- .Call("callNPSOL", objective, startVals, 
		constraints, matrices, parameters, 
		algebras, data, options, state, PACKAGE = "OpenMx")
	backendStop <- Sys.time()
	model <- omxReplaceModels(model, independents)
	model <- updateModelMatrices(model, flatModel, output$matrices)
	model <- updateModelAlgebras(model, flatModel, output$algebras)
	model <- undoSquareBracketLabels(model)
	model@output <- processOptimizerOutput(flatModel, names(matrices),
		names(algebras), names(parameters), output)
	frontendStop2 <- Sys.time()
	model@output <- calculateTiming(model@output, frontendStart,
		frontendStop1, backendStop, frontendStop2)	
	return(model)
}

processOptimizerOutput <- function(flatModel, matrixNames, 
		algebraNames, parameterNames, output) {
	output$mxVersion <- mxVersion()
	if (length(output$estimate) == length(parameterNames)) {
		names(output$estimate) <- parameterNames
	}
	if (length(output$gradient) == length(parameterNames)) {
		names(output$gradient) <- parameterNames
	}
	output$hessian <- t(output$hessianCholesky) %*% output$hessianCholesky
	if (nrow(output$hessian) == length(parameterNames) &&
		ncol(output$hessian) == length(parameterNames)) {
		dimnames(output$hessian) <- list(parameterNames, parameterNames)
	}
	if (length(output$matrices) == length(matrixNames)) {
		names(output$matrices) <- matrixNames
	}
	if (length(output$algebras) == length(algebraNames)) {
		names(output$algebras) <- algebraNames
	}
    if (output$status[[1]] > 0) {
    	npsolWarnings(flatModel@name, output$status[[1]])
    } else if (output$status[[1]] < 0) {
        stop(paste("The job for model", omxQuotes(flatModel@name),
            "exited abnormally with the error message:",
            output$status[[3]]), call. = FALSE)
    }
	return(output)
}

calculateTiming <- function(output, frontendStart,
	frontendStop1, backendStop, frontendStop2) {
	output$frontendTime <- (frontendStop1 - frontendStart) +
		(frontendStop2 - backendStop)
	output$backendTime <- (backendStop - frontendStop1)
	return(output)
}

npsolMessages <- list('1' = paste('The final iterate satisfies',
		'the optimality conditions to the accuracy requested,',
		'but the sequence of iterates has not yet converged.',
		'NPSOL was terminated because no further improvement',
		'could be made in the merit function (Mx status GREEN).'),
		'2' = paste('The linear constraints and bounds could not be satisfied.',
		'The problem has no feasible solution.'),
		'3' = paste('The nonlinear constraints and bonuds could not be satisfied.',
		'The problem may have no feasible solution.'),
		'4' = 'The Major iteration limit was reached (Mx status BLUE).',
		'6' = paste('model does not satisfy the first-order optimality conditions',
		'to the required accuracy, and no improved point for the',
		'merit function could be found during the final linesearch (Mx status RED)'),
		'7' = paste('The function derivates returned by funcon or funobj',
		'appear to be incorrect.'),
		'9' = 'An input parameter was invalid')

npsolWarnings <- function(name, status) {
	message <- npsolMessages[[as.character(status)]]
	if(!is.null(message)) {
		warning(paste("In model", omxQuotes(name), 
			"NPSOL returned a non-zero status code", 
			paste(status, '.', sep = ''), message), call. = FALSE)
	}
}