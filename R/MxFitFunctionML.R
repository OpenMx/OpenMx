#
#   Copyright 2007-2014 The OpenMx Project
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
	contains = "MxBaseFitFunction")

setMethod("initialize", "MxFitFunctionML",
	function(.Object, vector, name = 'fitfunction') {
		.Object@name <- name
		.Object@vector <- vector
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
		}

		return(flatModel)
})

setMethod("genericFitFunConvert", "MxFitFunctionML", 
	function(.Object, flatModel, model, labelsData, defVars, dependencies) {
		name <- .Object@name
		modelname <- imxReverseIdentifier(model, .Object@name)[[1]]
		expectName <- paste(modelname, "expectation", sep=".")
		if (expectName %in% names(flatModel@expectations)) {
			expectIndex <- imxLocateIndex(flatModel, expectName, name)
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
	function(.Object, model) {
		modelName <- model@name
		datasource <- model$data
		if (is.null(datasource)) {
			stop(paste("Model", omxQuotes(modelName), "does not contain any data"))
		}

		expectation <- model@expectation
		if (is(expectation, "MxExpectationBA81")) {
			return(generateIFAReferenceModels(model))
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
				obsdata <- obsdata[selVars, selVars]
			} else { obsdata <- obsdata[,selVars] }
		} else {
			message(paste("The model", omxQuotes(modelName), "has not been run. So reference models",
				"of all the variables in the data will be made.  For reference models",
				"of only the variables used in the model, provide the model after it has been run."))
		}
		if(.hasSlot(model@expectation, 'definitionVars')){
			if(length(model@expectation@definitionVars) > 0){
				message(paste("The model", omxQuotes(modelName), "has definition variables.",
				"The saturated and independence models generated may be incorrect.",
				"Commence shooting yourself in the foot."))
			}
		}
		
		generateNormalReferenceModels(modelName, obsdata, datatype, any(!is.na(datasource@means)), datanobs)
	})

mxFitFunctionML <- function(vector = FALSE) {
	if (length(vector) > 1 || typeof(vector) != "logical") {
		stop("'vector' argument is not a logical value")
	}
	return(new("MxFitFunctionML", vector))
}

displayMxFitFunctionML <- function(fitfunction) {
	cat("MxFitFunctionML", omxQuotes(fitfunction@name), '\n')
	cat("$vector :", fitfunction@vector, '\n')
	print(fitfunction@result)
	invisible(fitfunction)
}


setMethod("print", "MxFitFunctionML", function(x, ...) { 
	displayMxFitFunctionML(x) 
})

setMethod("show", "MxFitFunctionML", function(object) { 
	displayMxFitFunctionML(object) 
})
