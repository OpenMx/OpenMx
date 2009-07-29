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

mxEvaluate <- function(expression, model, show = FALSE) {
	formula <- match.call()$expression
	modelVariable <- match.call()$model
	labelsData <- omxGenerateLabels(model)
	formula <- evaluateTranslation(formula, model, modelVariable, labelsData, FALSE)
	if (show) {
		showFormula <- evaluateTranslation(formula, model, modelVariable, labelsData, TRUE)
		cat(deparse(showFormula, width.cutoff = 500L), '\n')
	}
	return(eval(formula))
}

evaluateTranslation <- function(formula, model, modelVariable, labelsData, show) {
	len <- length(formula)
	if (len == 0) {
		stop("mxEvaluate has reached an invalid state")
	} else if (len == 1) {
		formula <- translateSymbol(formula, model, modelVariable, labelsData, show)
	} else {
		formula[-1] <- lapply(formula[-1], 
			evaluateTranslation, model, modelVariable, labelsData, show)
	}
	return(formula)
}

translateSymbol <- function(symbol, model, modelVariable, labelsData, show) {
	key <- deparse(symbol)
	index <- match(key, dimnames(labelsData)[[1]])
	if (!show) modelVariable <- as.symbol('model')
	if (!is.na(index)) {
		labelModel <- labelsData[[index,"model"]]
		labelMatrix <- labelsData[[index,"matrix"]]
		labelRow <- labelsData[[index,"row"]]
		labelCol <- labelsData[[index,"col"]]
		return(substitute(modelName[[x]]@values[[y,z]],
			list(modelName = modelVariable,
				x = omxIdentifier(labelModel, labelMatrix),
				y = labelRow, z = labelCol)))
	}
	lookup <- model[[key]]
	if (is.null(lookup)) {
		return(symbol)
	} else if (is(lookup, "MxMatrix")) {
		return(substitute(modelName[[x]]@values,
			list(modelName = modelVariable, x = key)))
	} else if (is(lookup, "MxAlgebra")) {
		return(substitute(modelName[[x]]@result,
			list(modelName = modelVariable, x = key)))
	} else if (is(lookup, "MxObjective")) {
		return(substitute(modelName[[x]]@result,
			list(modelName = modelVariable, x = key)))
	} else {
		stop(paste("Cannot evaluate the object",
			omxQuotes(key), "in the model",
			omxQuotes(model@name)))
	}
}

omxGenerateLabels <- function(model) {
	return(generateLabelsHelper(model, data.frame()))
}

generateLabelsHelper <- function(model, labelsData) {
	if (length(model@matrices) > 0) {
		for(i in 1:length(model@matrices)) {
			labelsData <- generateLabelsMatrix(model@name, model@matrices[[i]], labelsData)
		}
	}
	if (length(model@submodels) > 0) {
		for(i in 1:length(model@submodels)) {
			labelsData <- generateLabelsHelper(model@submodels[[i]], labelsData)
		}
	}
	return(labelsData)
}

generateLabelsMatrix <- function(modelName, matrix, labelsData) {
	labels <- matrix@labels
	select <- labels[!is.na(labels)]
	rows <- row(labels)[!is.na(labels)]
	cols <- col(labels)[!is.na(labels)]
	if (length(select) > 0) {
		for(i in 1:length(select)) {
			if(!omxIsDefinitionVariable(select[[i]])) {
				labelsData[select[[i]], "model"] <- modelName
				labelsData[select[[i]], "matrix"] <- matrix@name
				labelsData[select[[i]], "row"] <- rows[[i]]
				labelsData[select[[i]], "col"] <- cols[[i]]
			}
		}
	}
	return(labelsData)
}
