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

# The mxEval funcion is NOT reentrant.  This is because the 
# convertSquareBracketLabels() transformation can only occur once.
# Recursive calls should use the helper function computeSymbol().
# There would also be a performance penalty associated with a 
# recursive call to mxEval().
mxEval <- function(expression, model, compute = FALSE, show = FALSE) {
	model <- convertSquareBracketLabels(model)
	inputExpression <- match.call()$expression
	labelsData <- omxGenerateLabels(model)
    if (compute) {
        namespace <- omxGenerateNamespace(model)
        flatModel <- omxFlattenModel(model, namespace)
        flatExpression <- namespaceConvertFormula(inputExpression, model@name, namespace)
    	formula <- computeTranslation(flatExpression, flatModel, labelsData)
    } else {
    	formula <- evaluateTranslation(inputExpression, model, labelsData)
    }
	if (show) {
    	modelVariable <- match.call()$model
        if (compute) {
            showFormula <- showEvaluation(flatExpression, flatModel, modelVariable, labelsData)
        } else {
    		showFormula <- showTranslation(inputExpression, model, modelVariable, labelsData)
        }
		cat(deparse(showFormula, width.cutoff = 500L), '\n')
	}
	return(eval(formula))
}

evaluateTranslation <- function(formula, model, labelsData) {
	len <- length(formula)
	if (len == 0) {
		stop("mxEval has reached an invalid state")
	} else if (len == 1) {
		formula <- translateSymbol(formula, model, labelsData)
	} else if (len == 4 && identical(as.character(formula[1]), '[')) {
		formula[-1] <- lapply(formula[-1], 
			evaluateTranslation, model, labelsData)
		# The order of the if/else-if matters
		# as.matrix(x) should be invoked on A[,y] and A[,]
		# t(as.matrix(x)) should be invoked on A[x,]
		if (identical(as.character(formula[3]), '')) {
			formula <- substitute(as.matrix(x), list(x = formula))
		} else if (identical(as.character(formula[4]), '')) {
			formula <- substitute(t(as.matrix(x)), list(x = formula))
		}
	} else {
		formula[-1] <- lapply(formula[-1], 
			evaluateTranslation, model, labelsData)
	}
	return(formula)
}

translateSymbol <- function(symbol, model, labelsData) {
	key <- deparse(symbol)
	index <- match(key, dimnames(labelsData)[[1]])
	if (!is.na(index)) {
		labelModel <- labelsData[[index,"model"]]
		labelMatrix <- labelsData[[index,"matrix"]]
		labelRow <- labelsData[[index,"row"]]
		labelCol <- labelsData[[index,"col"]]
		return(substitute(model[[x]]@values[[y,z]],
			list(x = omxIdentifier(labelModel, labelMatrix),
				y = labelRow, z = labelCol)))
	}
	lookup <- model[[key]]
	if (is.null(lookup)) {
		return(symbol)
	} else if (is(lookup, "MxMatrix")) {
		return(substitute(model[[x]]@values, list(x = key)))
	} else if (is(lookup, "MxAlgebra")) {
		return(substitute(model[[x]]@result, list(x = key)))
	} else if (is(lookup, "MxObjective")) {
		return(substitute(model[[x]]@result, list(x = key)))
	} else {
		stop(paste("Cannot evaluate the object",
			omxQuotes(key), "in the model",
			omxQuotes(model@name)))
	}
}

computeTranslation <- function(formula, model, labelsData) {
	len <- length(formula)
	if (len == 0) {
		stop("mxEval has reached an invalid state")
	} else if (len == 1) {
		formula <- computeSymbol(formula, model, labelsData)
	} else if (len == 4 && identical(as.character(formula[1]), '[')) {
		formula[-1] <- lapply(formula[-1], 
			computeTranslation, model, labelsData)
		# The order of the if/else-if matters
		# as.matrix(x) should be invoked on A[,y] and A[,]
		# t(as.matrix(x)) should be invoked on A[x,]
		if (identical(as.character(formula[3]), '')) {
			formula <- substitute(as.matrix(x), list(x = formula))
		} else if (identical(as.character(formula[4]), '')) {
			formula <- substitute(t(as.matrix(x)), list(x = formula))
		}
	} else {
		formula[-1] <- lapply(formula[-1], 
			computeTranslation, model, labelsData)
	}
	return(formula)
}

captureComputeTranslation <- function(expression, model, labelsData) {
	formula <- match.call()$expression
	return(computeTranslation(formula, model, labelsData))
}

computeSymbol <- function(symbol, model, labelsData) {
	key <- deparse(symbol)
	index <- match(key, dimnames(labelsData)[[1]])
	if (!is.na(index)) {
		labelModel <- labelsData[[index,"model"]]
		labelMatrix <- labelsData[[index,"matrix"]]
		labelRow <- labelsData[[index,"row"]]
		labelCol <- labelsData[[index,"col"]]
		return(substitute(model[[x]]@values[[y,z]],
			list(x = omxIdentifier(labelModel, labelMatrix),
				y = labelRow, z = labelCol)))
	}
	lookup <- model[[key]]
	if (omxIsDefinitionVariable(key)) {
		return(definitionStartingValue(key, model))
	} else if (is.null(lookup)) {
		return(symbol)
	} else if (is(lookup, "MxMatrix")) {
		return(substitute(computeMatrixHelper(x, flatModel, labelsData), list(x = lookup)))
	} else if (is(lookup, "MxAlgebra")) {
		return(substitute(omxDimnames(eval(captureComputeTranslation(x, flatModel, labelsData)), y),
			list(x = lookup@formula, y = lookup@.dimnames)))
	} else if (is(lookup, "MxObjective")) {
        if (length(lookup@result) == 0) {
            return(omxObjInitialMatrix(lookup, model))
        } else {
    		return(substitute(model[[x]]@result, list(x = key)))
        }
	} else {
		stop(paste("Cannot evaluate the object",
			omxQuotes(key), "in the model",
			omxQuotes(model@name)))
	}
}

omxComputeMatrix <- function(matrix, model) {
	return(eval(substitute(mxEval(x, model, TRUE), list(x = as.symbol(matrix@name)))))
}

computeMatrixHelper <- function(matrix, flatModel, labelsData) {
	labels <- matrix@labels
	select <- !apply(labels, c(1,2), is.na) & apply(labels, c(1,2), hasSquareBrackets)
	subs <- labels[select]
	rows <- row(labels)[select]
	cols <- col(labels)[select]
	if (length(subs) == 0) {
		return(matrix@values)
	}
	for (i in 1:length(subs)) {
		algname <- subs[[i]]
		result <- tryCatch(eval(computeSymbol(as.symbol(algname), flatModel, labelsData)),
			error = function(x) {
				algebra <- flatModel[[algname]]
				stop(paste("The label", 
					omxQuotes(simplifyName(deparse(algebra@formula, width.cutoff=500L), flatModel@name)),
					"of matrix", omxQuotes(simplifyName(matrix@name, flatModel@name)),
					"in model", omxQuotes(flatModel@name), 
					"generated the error message:",
					x$message), call. = FALSE)
		})
		result <- as.matrix(result)
		if (nrow(result) != 1 || ncol(result) != 1) {
			algebra <- flatModel[[algname]]
			stop(paste("The label", 
				omxQuotes(simplifyName(deparse(algebra@formula, width.cutoff=500L), flatModel@name)),
				"of matrix", omxQuotes(simplifyName(matrix@name, flatModel@name)),
				"in model", omxQuotes(flatModel@name), 
				"does not evaluate to a (1 x 1) matrix."), call. = FALSE)
		}
		matrix@values[rows[[i]], cols[[i]]] <- result[1,1]
	}
	return(matrix@values)
}

omxDimnames <- function(value, names) {
	value <- as.matrix(value)
	dimnames(value) <- names
	return(value)
}

showTranslation <- function(formula, model, modelVariable, labelsData) {
	len <- length(formula)
	if (len == 0) {
		stop("mxEval has reached an invalid state")
	} else if (len == 1) {
		formula <- showSymbol(formula, model, modelVariable, labelsData)
	} else if (len == 4 && identical(as.character(formula[1]), '[')) {
		formula[-1] <- lapply(formula[-1], 
			showTranslation, model, modelVariable, labelsData)
		# The order of the if/else-if matters
		# as.matrix(x) should be invoked on A[,y] and A[,]
		# t(as.matrix(x)) should be invoked on A[x,]
		if (identical(as.character(formula[3]), '')) {
			formula <- substitute(as.matrix(x), list(x = formula))
		} else if (identical(as.character(formula[4]), '')) {
			formula <- substitute(t(as.matrix(x)), list(x = formula))
		}
	} else {
		formula[-1] <- lapply(formula[-1], 
			showTranslation, model, modelVariable, labelsData)
	}
	return(formula)
}

showSymbol <- function(symbol, model, modelVariable, labelsData) {
	key <- deparse(symbol)
	index <- match(key, dimnames(labelsData)[[1]])
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

showEvaluation <- function(formula, model, modelVariable, labelsData) {
	len <- length(formula)
	if (len == 0) {
		stop("mxEval has reached an invalid state")
	} else if (len == 1) {
		formula <- showEvaluationSymbol(formula, model, modelVariable, labelsData)
	} else if (len == 4 && identical(as.character(formula[1]), '[')) {
		formula[-1] <- lapply(formula[-1], 
			showEvaluation, model, modelVariable, labelsData)
		# The order of the if/else-if matters
		# as.matrix(x) should be invoked on A[,y] and A[,]
		# t(as.matrix(x)) should be invoked on A[x,]
		if (identical(as.character(formula[3]), '')) {
			formula <- substitute(as.matrix(x), list(x = formula))
		} else if (identical(as.character(formula[4]), '')) {
			formula <- substitute(t(as.matrix(x)), list(x = formula))
		}
	} else {
		formula[-1] <- lapply(formula[-1], 
			showEvaluation, model, modelVariable, labelsData)
	}
	return(formula)
}

showEvaluationSymbol <- function(symbol, model, modelVariable, labelsData) {
	key <- deparse(symbol)
	index <- match(key, dimnames(labelsData)[[1]])
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
		labels <- lookup@labels
		select <- !apply(labels, c(1,2), is.na) & apply(labels, c(1,2), hasSquareBrackets)
	    subs <- labels[select]
	    if (length(subs) == 0) {
			return(substitute(modelName[[x]]@values, 
				list(modelName = modelVariable, x = key)))
		} else {
			return(substitute(omxComputeMatrix(modelName[[x]], model),
				list(modelName = modelVariable, x = key)))
		}
	} else if (is(lookup, "MxAlgebra")) {
		return(substitute(mxEval(x, modelName, TRUE), list(x = lookup@formula,
            modelName = modelVariable)))
	} else if (is(lookup, "MxObjective")) {
        if (length(lookup@result) == 0) {
            return(omxObjInitialMatrix(lookup, model))
        } else {
    		return(substitute(modelName[[x]]@result, 
                list(x = key, modelName = modelVariable)))
        }
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
			if(!omxIsDefinitionVariable(select[[i]]) && !hasSquareBrackets(select[[i]])) {
				labelsData[select[[i]], "model"] <- modelName
				labelsData[select[[i]], "matrix"] <- matrix@name
				labelsData[select[[i]], "row"] <- rows[[i]]
				labelsData[select[[i]], "col"] <- cols[[i]]
			}
		}
	}
	return(labelsData)
}
