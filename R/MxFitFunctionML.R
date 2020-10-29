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


setClass(Class = "MxFitFunctionML",
	 contains = "MxBaseFitFunction",
	 representation = representation(
	     fellner = "logical",
	     verbose = "integer",
	     profileOut="MxOptionalChar",
	     rowwiseParallel="logical",
	     jointConditionOn="character",
	     components="MxCharOrNumber"),
	 )

setMethod("initialize", "MxFitFunctionML",
	  function(.Object, vector, rowDiagnostics, fellner, verbose, profileOut,
		   rowwiseParallel, jointConditionOn, name = 'fitfunction') {
		.Object@name <- name
		.Object@vector <- vector
		.Object@rowDiagnostics <- rowDiagnostics
		.Object@fellner <- fellner
		.Object@verbose <- verbose
		.Object@profileOut <- profileOut
		.Object@rowwiseParallel <- rowwiseParallel
		.Object@jointConditionOn <- jointConditionOn
		return(.Object)
	}
)

setMethod("qualifyNames", signature("MxFitFunctionML"), 
	function(.Object, modelname, namespace) {
		.Object@name <- imxIdentifier(modelname, .Object@name)
		return(.Object)
})

setMethod("genericFitConvertEntities", "MxFitFunctionML",
	function(.Object, flatModel, namespace, labelsData) {

		name <- .Object@name
		modelname <- imxReverseIdentifier(flatModel, .Object@name)[[1]]
		expectName <- paste(modelname, "expectation", sep=".")

		expectation <- flatModel@expectations[[expectName]]
		dataname <- expectation@data		

		if (flatModel@datasets[[dataname]]@type != 'raw') {
			if (.Object@vector) {
				modelname <- getModelName(.Object)
				msg <- paste("The ML fit function",
					"in model", omxQuotes(modelname), "has specified",
					"'vector' = TRUE, but the observed data is not raw data")
				stop(msg, call.=FALSE)
			}
			if (.Object@rowDiagnostics) {
				modelname <- getModelName(.Object)
				msg <- paste("The ML fit function",
					"in model", omxQuotes(modelname), "has specified",
					"'rowDiagnostics' = TRUE, but the observed data is not raw data")
				stop(msg, call.=FALSE)
			}
		}

		return(flatModel)
})

setMethod("genericFitFunConvert", "MxFitFunctionML", 
	function(.Object, flatModel, model, labelsData, dependencies) {
		name <- .Object@name
		modelname <- imxReverseIdentifier(model, .Object@name)[[1]]
		expectName <- paste(modelname, "expectation", sep=".")
		if (expectName %in% names(flatModel@expectations)) {
			expectIndex <- imxLocateIndex(flatModel, expectName, name)

			ex <- flatModel@expectations[[1L + expectIndex]]
			if (is(ex, "MxExpectationHiddenMarkov") || is(ex, "MxExpectationMixture")) {
				.Object@components <-
					sapply(paste(ex@components, "fitfunction", sep="."),
					       function(ff) imxLocateIndex(flatModel, ff, name),
					       USE.NAMES = FALSE)
			}
		} else {
			expectIndex <- as.integer(NA)
		}
		.Object@expectation <- expectIndex
		return(.Object)
})


setMethod("genericFitInitialMatrix", "MxFitFunctionML",
	function(.Object, flatModel) {
		flatFitFunction <- flatModel@fitfunctions[[.Object@name]]
		if (flatFitFunction@vector == FALSE) {
			return(matrix(as.double(NA), 1, 1))
		} else {
			modelname <- imxReverseIdentifier(flatModel, flatFitFunction@name)[[1]]
			expectationName <- paste(modelname, "expectation", sep = ".")
			expectation <- flatModel@expectations[[expectationName]]
			if (is.null(expectation)) {
				msg <- paste("The ML fit function has vector = TRUE",
				"and a missing expectation in the model",
				omxQuotes(modelname))
				stop(msg, call.=FALSE)
			}
			if (is.na(expectation@data)) {
				msg <- paste("The ML fit function has vector = TRUE",
				"and an expectation function with no data in the model",
				omxQuotes(modelname))
				stop(msg, call.=FALSE)
			}
			mxDataObject <- flatModel@datasets[[expectation@data]]
			if (mxDataObject@type != 'raw') {
				msg <- paste("The dataset associated with the ML expectation function", 
					"in model", omxQuotes(modelname), "is not raw data.")
				stop(msg, call.=FALSE)
			}
			rows <- nrow(mxDataObject@observed)
			return(matrix(as.double(NA), rows, 1))
		}
})

setMethod("generateReferenceModels", "MxFitFunctionML",
	function(.Object, model, distribution, equateThresholds) {
		modelName <- model@name
		datasource <- model$data
		if (is.null(datasource)) {
			stop(paste("Model", omxQuotes(modelName), "does not contain any data"))
		}

		expectation <- model@expectation
		if (is(expectation, "MxExpectationBA81")) {
			return(generateIFAReferenceModels(model, distribution))
		}
		if(is(expectation, "MxExpectationGREML")){
			stop("Reference models for GREML expectation are not implemented")
		}
		# assume it's multivariate Normal

		datatype <- datasource@type
		obsdata <- datasource@observed
		datanobs <- datasource@numObs
		wasRun <- model@.wasRun
		if(wasRun) {
			if (is.null(model@expectation@.runDims)) stop("Not clear which data were used to fit model")
			selVars <- model@expectation@.runDims
			if(nrow(obsdata) == ncol(obsdata)){
				if(!single.na(model@expectation@.runDims)) { obsdata <- obsdata[selVars, selVars] }
				#variable subsets are not run for covariance data
				#consequently, selVars are only used when runDims are provided.
			} else { obsdata <- obsdata[,selVars, drop=FALSE] }
		} else {
			message(paste("The model", omxQuotes(modelName), "has not been run. So reference models",
				"of all the variables in the data will be made.  For reference models",
				"of only the variables used in the model, provide the model after it has been run."))
		}
		
		generateNormalReferenceModels(modelName, obsdata, datatype, any(!is.na(datasource@means)),
			datanobs, datasource@means, distribution=distribution,
			equateThresholds)
	})

mxFitFunctionML <- function(vector = FALSE, rowDiagnostics=FALSE, ..., fellner=as.logical(NA),
			    verbose=0L, profileOut=c(), rowwiseParallel=as.logical(NA),
			    jointConditionOn=c('auto', 'ordinal', 'continuous')) {
	prohibitDotdotdot(list(...))
	if (length(vector) > 1 || typeof(vector) != "logical") {
		stop("'vector' argument is not a logical value")
	}
	if (length(rowDiagnostics) > 1 || typeof(rowDiagnostics) != "logical") {
		stop("'rowDiagnostics' argument is not a logical value")
	}
	if (length(fellner) > 1) {
		stop("'fellner' argument must be one thing")
	}
	if (!is.na(fellner) && fellner && (vector || rowDiagnostics)) {
		stop("'fellner' cannot be combined with 'vector' or 'rowDiagnostics'")
	}
	jointConditionOn <- match.arg(jointConditionOn)
	return(new("MxFitFunctionML", vector, rowDiagnostics, fellner,
		   as.integer(verbose), as.character(profileOut), rowwiseParallel,
		   jointConditionOn))
}

displayMxFitFunctionML <- function(fitfunction) {
	cat("MxFitFunctionML", omxQuotes(fitfunction@name), '\n')
	cat("$vector :", fitfunction@vector, '\n')
	cat("$rowDiagnostics :", fitfunction@rowDiagnostics, '\n')
	cat("$fellner :", fitfunction@fellner, '\n')
	cat("$verbose :", fitfunction@verbose, '\n')
	cat("$rowwiseParallel :", fitfunction@rowwiseParallel, '\n')
	cat("$jointConditionOn :", fitfunction@jointConditionOn, '\n')
	print(fitfunction@result)
	invisible(fitfunction)
}


setMethod("print", "MxFitFunctionML", function(x, ...) { 
	displayMxFitFunctionML(x) 
})

setMethod("show", "MxFitFunctionML", function(object) { 
	displayMxFitFunctionML(object) 
})
