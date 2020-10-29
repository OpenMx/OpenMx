#
#   Copyright 2007-2019 by the individuals mentioned in the source code history
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
		algebraNames, parameterNames, constraintNames, computeSeq, 
		output) {
	output$mxVersion <- mxVersion(verbose=FALSE)
	if (length(output$estimate) == length(parameterNames)) {
		names(output$estimate) <- parameterNames
	}
	for (deriv in c("calculatedHessian")) {
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
	output <- nameGenericConstraintOutput(parameterNames, constraintNames, output)
	GDsteps <- computeSeq$steps[sapply(computeSeq$steps,function(x){is(x,"MxComputeGradientDescent")})]
	if(length(GDsteps)){
		lastGDstep <- GDsteps[[length(GDsteps)]] #<--If a custom compute plan has more than one GD step, use the last one
		output <- nameGDOptimizerConstraintOutput(parameterNames, constraintNames, lastGDstep, output)
	}
	return(output)
}


nameGDOptimizerConstraintOutput <- function(paramNames, constraintNames, GDstep, output){
	if(GDstep@engine=="NPSOL"){
		
		#Initialize variables:
		cfvNames <- NULL
		lmNames <- NULL
		constraintRows <- NULL
		constraintCols <- NULL
		
		#Filter extraneous elements and generate vectors of names:
		if(length(paramNames)){lmNames <- paste(paramNames,"bound",sep=".")}
		if(length(constraintNames) && length(GDstep@output$constraintRows) && length(GDstep@output$constraintCols)){
			emptyConstraints <- (GDstep@output$constraintRows==0 | GDstep@output$constraintCols==0)
			#Assuming that "empty" constraints have no function values...
			constraintNames <- constraintNames[!emptyConstraints]
			constraintRows <- GDstep@output$constraintRows[!emptyConstraints]
			constraintCols <- GDstep@output$constraintCols[!emptyConstraints]
		}
		if(length(constraintNames) && length(constraintRows) && length(constraintCols)){
			for(i in 1:length(constraintNames)){
				for(co in 1:constraintCols[i]){
					for(ro in 1:constraintRows[i]){
						cfvNames <- c(cfvNames, paste(constraintNames[i],"[",ro,",",co,"]",sep=""))
					}
				}
			}
			lmNames <- c(lmNames, cfvNames)
		}
		
		#Assign names to components, and components to 'output' list:
		if(is.null(output$constraintFunctionValues)){
			if(length(GDstep@output$constraintFunctionValues)){
				if(length(GDstep@output$constraintFunctionValues)==length(cfvNames)){
					names(GDstep@output$constraintFunctionValues) <- cfvNames
				}
				output$constraintFunctionValues <- GDstep@output$constraintFunctionValues
			}
		}
		if(is.null(output$constraintJacobian)){
			if(length(GDstep@output$constraintJacobian)){
				if(ncol(GDstep@output$constraintJacobian) && ncol(GDstep@output$constraintJacobian)==length(paramNames)){
					colnames(GDstep@output$constraintJacobian) <- paramNames
				}
				if(nrow(GDstep@output$constraintJacobian) && nrow(GDstep@output$constraintJacobian)==length(cfvNames)){
					rownames(GDstep@output$constraintJacobian) <- cfvNames
				}
				output$constraintJacobian <- GDstep@output$constraintJacobian
			}
		}
		if(length(GDstep@output$LagrangeMultipliers)){
			if(length(GDstep@output$LagrangeMultipliers)==length(lmNames)){
				names(GDstep@output$LagrangeMultipliers) <- lmNames
			}
			output$LagrangeMultipliers <- GDstep@output$LagrangeMultipliers
		}
		if(length(GDstep@output$istate)){
			if(length(GDstep@output$istate)==length(lmNames)){
				names(GDstep@output$istate) <- lmNames
			}
			output$istate <- GDstep@output$istate
		}
	}
	else if(GDstep@engine=="CSOLNP" || GDstep@engine=="SLSQP"){
		
		#Initialize variables:
		cfvNames <- NULL
		constraintRows <- NULL
		constraintCols <- NULL
		
		#Filter extraneous elements and generate vectors of names:
		if(length(constraintNames) && length(GDstep@output$constraintRows) && length(GDstep@output$constraintCols)){
			emptyConstraints <- (GDstep@output$constraintRows==0 | GDstep@output$constraintCols==0)
			#Assuming that "empty" constraints have no function values...
			constraintNames <- constraintNames[!emptyConstraints]
			constraintRows <- GDstep@output$constraintRows[!emptyConstraints]
			constraintCols <- GDstep@output$constraintCols[!emptyConstraints]
		}
		if(length(constraintNames) && length(constraintRows) && length(constraintCols)){
			for(i in 1:length(constraintNames)){
				for(co in 1:constraintCols[i]){
					for(ro in 1:constraintRows[i]){
						cfvNames <- c(cfvNames, paste(constraintNames[i],"[",ro,",",co,"]",sep=""))
					}
				}
			}
		}
		
		#Assign names to components, and components to 'output' list:
		if(is.null(output$constraintFunctionValues)){
			if(length(GDstep@output$constraintFunctionValues)){
				if(length(GDstep@output$constraintFunctionValues)==length(cfvNames)){
					names(GDstep@output$constraintFunctionValues) <- cfvNames
				}
				output$constraintFunctionValues <- GDstep@output$constraintFunctionValues
			}
		}
		if(is.null(output$constraintJacobian)){
			if(length(GDstep@output$constraintJacobian)){
				if(ncol(GDstep@output$constraintJacobian) && ncol(GDstep@output$constraintJacobian)==length(paramNames)){
					colnames(GDstep@output$constraintJacobian) <- paramNames
				}
				if(nrow(GDstep@output$constraintJacobian) && nrow(GDstep@output$constraintJacobian)==length(cfvNames)){
					rownames(GDstep@output$constraintJacobian) <- cfvNames
				}
				output$constraintJacobian <- GDstep@output$constraintJacobian
			}
		}
		if(length(GDstep@output$LagrangeMultipliers)){
			if(length(GDstep@output$LagrangeMultipliers)==length(cfvNames)){
				names(GDstep@output$LagrangeMultipliers) <- cfvNames
			}
			output$LagrangeMultipliers <- GDstep@output$LagrangeMultipliers
		}
		if(length(GDstep@output$LagrHessian)){
			if(ncol(GDstep@output$LagrHessian) && ncol(GDstep@output$LagrHessian)==length(paramNames) && 
				 nrow(GDstep@output$LagrHessian)==length(paramNames)){
				dimnames(GDstep@output$LagrHessian) <- list(paramNames,paramNames)
			}
			output$LagrHessian <- GDstep@output$LagrHessian
		}
	}
	return(output)
}


nameGenericConstraintOutput <- function(paramNames, constraintNames, output){
	#Name columns of vcov, if it exists:
	if(!is.null(output$vcov)){ 
		 if(ncol(output$vcov) && nrow(output$vcov) && nrow(output$vcov)==ncol(output$vcov)){
			colnames(output$vcov) <- rownames(output$vcov)
	}}
	
	#Initialize variables:
	cfvNames <- NULL
	constraintRows <- NULL
	constraintCols <- NULL
	
	#Filter extraneous elements and generate vectors of names:
	if(length(constraintNames) && length(output$constraintRows) && length(output$constraintCols)){
		emptyConstraints <- (output$constraintRows==0 | output$constraintCols==0)
		#Assuming that "empty" constraints have no function values...
		constraintNames <- constraintNames[!emptyConstraints]
		constraintRows <- output$constraintRows[!emptyConstraints]
		constraintCols <- output$constraintCols[!emptyConstraints]
	}
	if(length(constraintNames) && length(constraintRows) && length(constraintCols)){
		for(i in 1:length(constraintNames)){
			for(co in 1:constraintCols[i]){
				for(ro in 1:constraintRows[i]){
					cfvNames <- c(cfvNames, paste(constraintNames[i],"[",ro,",",co,"]",sep=""))
				}
			}
		}
	}
	
	if(length(output$constraintFunctionValues) && length(output$constraintFunctionValues)==length(cfvNames)){
		names(output$constraintFunctionValues) <- cfvNames
	}
	if(length(output$constraintJacobian)){
		if(ncol(output$constraintJacobian) && ncol(output$constraintJacobian)==length(paramNames)){
			colnames(output$constraintJacobian) <- paramNames
		}
		if(nrow(output$constraintJacobian) && nrow(output$constraintJacobian)==length(cfvNames)){
			rownames(output$constraintJacobian) <- cfvNames
		}
	}
	
	#Remove clutter from output list:
	output$constraintNames <- NULL
	output$constraintCols <- NULL
	output$constraintRows <- NULL
	
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
		'5' = 'The Hessian at the solution does not appear to be convex. See ?mxCheckIdentification for possible diagnosis (Mx status RED).',
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
