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

mxEval <- function(expression, model, compute = FALSE, show = FALSE, defvar.row = 1) {
	if (missing(expression)) {
		stop("'expression' argument is mandatory in call to mxEval function")
	} else if (missing(model)) {
		stop("'model' argument is mandatory in call to mxEval function")
	} else if (!is.numeric(defvar.row) || length(defvar.row) != 1 || is.na(defvar.row)) {
		stop("'defvar.row' argument must be a single numeric value")
	}	
	expression <- match.call()$expression
	modelvariable <- match.call()$model
	labelsData <- imxGenerateLabels(model)
	env <- parent.frame()
	if (compute) {
		namespace <- imxGenerateNamespace(model)
		model <- imxFlattenModel(model, namespace)
		expression <- namespaceConvertFormula(expression, model@name, namespace)
	}
	if (show) {
		showResult <- evaluateExpression(expression, deparse(expression), model, 
			labelsData, env, compute, show = TRUE, outsideAlgebra = TRUE)
		showResult <- eval(substitute(substitute(x, list(.zzz = modelvariable)), list(x = showResult)))
		print(deparse(showResult), width.cutoff = 450L)
	}
	result <- evaluateExpression(expression, deparse(expression), model, 
		labelsData, env, compute, show = FALSE, outsideAlgebra = TRUE, defvar.row)
	if (is.vector(result)) {
		return(result)
	} else {
		return(eval(result, envir = env))
	}
}

translateErrorSymbol <- function(symbol, model) {
	charsymbol <- as.character(symbol)
	object <- model[[charsymbol]]
	if (!is.null(object) && is(object, "MxMatrix") && length(object@display) == 1) {
		return(as.symbol(object@display))
	} else {
		return(symbol)
	}
}

translateErrorFormula <- function(formula, model) {
	if(length(formula) == 1) {
		formula <- translateErrorSymbol(formula, model)
	} else {
		formula <- lapply(formula, translateErrorFormula, model)
		formula <- as.call(formula)
	}
	return(formula)
}

evaluateExpression <- function(formula, contextString, model, labelsData, env, compute, show, outsideAlgebra, defvar.row = 1) {
	len <- length(formula)
	if (len == 0) {
		stop("mxEval has reached an invalid state")
	} else if (len == 1) {
		if (!identical(as.character(formula), "")) {
			formula <- evaluateSymbol(formula, contextString, model, 
				labelsData, env, compute, show, outsideAlgebra, defvar.row)
		}
		return(formula)
	}
	originalFormula <- formula
	formula[-1] <- lapply(formula[-1], 
		evaluateExpression, contextString, model, labelsData,
			env, compute, show, outsideAlgebra, defvar.row)
	if (len == 4 && identical(as.character(formula[1]), '[')) {
		formula$drop <- FALSE
		if(is.matrix(formula[[3]]) && nrow(formula[[3]]) == 0 && ncol(formula[[3]]) == 0) {
			formula[[3]] <- TRUE
		}
		if(is.matrix(formula[[4]]) && nrow(formula[[4]]) == 0 && ncol(formula[[4]]) == 0) {
			formula[[4]] <- TRUE
		}
	}
	if (len == 3 && identical(as.character(formula[1]), '*')) {
		formula[[2]] <- substitute(as.matrix(x), list(x = formula[[2]]))
		formula[[3]] <- substitute(as.matrix(x), list(x = formula[[3]]))
	}
	if (len == 3 && identical(as.character(formula[1]), '^')) {
		formula[[2]] <- substitute(as.matrix(x), list(x = formula[[2]]))
		formula[[3]] <- substitute(as.matrix(x), list(x = formula[[3]]))
	}
	if (show) {
		return(formula)
	}
	result <- tryCatch(do.call(as.character(formula[[1]]), as.list(formula[-1]), envir = env),
			error = function(x) {
				originalFormula <- translateErrorFormula(originalFormula, model)
				formulaString <- deparse(originalFormula)
				if(identical(formulaString, contextString)) {
					msg <- paste("The following error occurred while",
						"evaluating the expression", omxQuotes(formulaString), 
						"in model", omxQuotes(model@name),
						":", x$message)
				} else {
					msg <- paste("The following error occurred while",
						"evaluating the subexpression",
						omxQuotes(formulaString), "during the evaluation of",
						omxQuotes(contextString), "in model", omxQuotes(model@name),
						":", x$message)
				}
				stop(msg, call. = FALSE)
		})
	return(result)
}

evaluateSymbol <- function(symbol, contextString, model, labelsData,
			env, compute, show, outsideAlgebra, defvar.row = 1) {
	key <- deparse(symbol)
	index <- match(key, dimnames(labelsData)[[1]])
	if (!is.na(index)) {
		targetmodel <- labelsData[[index,"model"]]
		matrix <- labelsData[[index,"matrix"]]
		fullname <- imxIdentifier(targetmodel, matrix)
		row <- labelsData[[index,"row"]]
		col <- labelsData[[index,"col"]]
		value <- model[[fullname]]@values[row,col]
		if (show) {
			return(substitute(.zzz[[x]]@values[y,z], list(x = fullname, y = row, z = col)))
		} else {
			return(value)
		}
	}
	lookup <- model[[key]]
	if (imxIsDefinitionVariable(key)) {
		return(definitionStartingValue(key, contextString, model, defvar.row))
	} else if (is.null(lookup)) {
		if (!show && !outsideAlgebra && exists(key, envir = env)) {
			return(as.matrix(get(key, envir = env)))
		} else {
			return(symbol)
		}
	} else if (is(lookup, "MxMatrix")) {
		return(evaluateMatrix(lookup, model, labelsData, env, show, compute, defvar.row))
	} else if (is(lookup, "MxAlgebra")) {
		return(evaluateAlgebra(lookup, model, labelsData, env, show, compute, defvar.row))
	} else if (is(lookup, "MxData")) {
		if (show) {
			return(substitute(.zzz[[x]]@observed, list(x = key)))
		} else {
			return(lookup@observed)
		}
	} else if (is(lookup, "MxObjective")) {
        if (length(lookup@result) == 0) {
            return(genericObjInitialMatrix(lookup, model))
        } else if (show) {
			return(substitute(.zzz[[x]]@result, list(x = key)))
		} else {
    		return(lookup@result)
        }
	} else {
		stop(paste("Cannot evaluate the object",
			omxQuotes(key), "in the model",
			omxQuotes(model@name)))
	}
}

evaluateMatrix <- function(matrix, model, labelsData, env, show, compute, defvar.row) {
	if (show) {
		return(substitute(.zzz[[x]]@values, list(x = matrix@name)))
	} else if (compute) {
		result <- computeMatrix(matrix, model, labelsData, defvar.row, env)
	} else {
		result <- matrix@values
	}
	result <- assignDimnames(matrix, result)
	return(result)
}

computeMatrix <- function(matrix, model, labelsData, defvar.row, env) {
	values <- populateDefVarMatrix(matrix, model, defvar.row)
	labels <- matrix@labels
	select <- !apply(labels, c(1,2), is.na) & apply(labels, c(1,2), hasSquareBrackets)
	if (all(!select)) {
		return(values)
	}
	subs <- labels[select]
	rows <- row(labels)[select]
	cols <- col(labels)[select]
	for (i in 1:length(subs)) {
		substitution <- subs[[i]]
		row <- rows[[i]]
		col <- cols[[i]]
		components <- splitSubstitution(substitution)
		expression <- parse(text = paste('`', components[[1]], '`[',
			as.integer(components[[2]]), ',', as.integer(components[[3]]), ']', sep=''))[[1]]
		contextString <- paste("label at row ", row, " and column ", col, " of matrix ", 
			omxQuotes(simplifyName(matrix@name, model@name)), sep = '')
		result <- evaluateExpression(expression, contextString, model, labelsData, env, 
			compute=TRUE, show=FALSE, outsideAlgebra=FALSE)
		result <- as.matrix(result)
		if (nrow(result) != 1 || ncol(result) != 1) {
			stop(paste("The label", 
				omxQuotes(substitution),
				"of matrix", omxQuotes(simplifyName(matrix@name, model@name)),
				"in model", omxQuotes(model@name), 
				"does not evaluate to a (1 x 1) matrix."), call. = FALSE)
		}
		values[row, col] <- result[1,1]
	}
	return(values)
}

evaluateAlgebra <- function(algebra, model, labelsData, env, show, compute, defvar.row) {
	if (compute && show) {
		return(computeAlgebra(algebra, model, labelsData, show = TRUE, defvar.row, env))
	} else if (compute) {
		result <- computeAlgebra(algebra, model, labelsData, show = FALSE, defvar.row, env)
		result <- as.matrix(result)
	} else if (show) {
		return(substitute(.zzz[[x]]@result, list(x = algebra@name)))
	} else {
		result <- algebra@result
	}
	if (!is.null(dimnames(algebra))) {
		result <- assignDimnames(algebra, result)
	}
	return(result)
}

computeAlgebra <- function(algebra, model, labelsData, show, defvar.row, env) {
	contextString <- simplifyName(algebra@name, model@name)
	result <- evaluateExpression(algebra@formula, contextString, model, labelsData, env, 
		compute = TRUE, show, outsideAlgebra = FALSE, defvar.row)
	return(result)
}

evaluateMxObject <- function(objname, flatModel, labelsData) {
	return(eval(substitute(evaluateSymbol(x, objname, flatModel, 
			labelsData, globalenv(), compute = TRUE, 
			show = FALSE, outsideAlgebra = FALSE, defvar.row = 1),
			list(x = quote(as.symbol(objname))))))
}

evaluateAlgebraWithContext <- function(algebra, context, flatModel, labelsData) {
	return(evaluateExpression(algebra@formula, context, flatModel, labelsData, globalenv(), 
		compute = TRUE, show = FALSE, outsideAlgebra = FALSE))
}

getEntityType <- function(object) {
	if(is(object, "MxMatrix")) {
		entity <- "MxMatrix"
	} else if(is(object, "MxAlgebra")) {
		entity <- "MxAlgebra"
	} else if(is(object, "MxObjective")) {
		entity <- "MxObjective"
	} else {
		entity <- "entity"
	}
	return(entity)
}

assignDimnames <- function(object, values) {
	newnames <- dimnames(object)
	if (!is.null(newnames)) {
		if (length(newnames) != 2) {
			entity <- getEntityType(object)
			msg <- paste("The 'dimnames' argument to", entity,
				omxQuotes(object@name), "must have a length of 2")
			stop(msg, call. = FALSE)
		}
		if (!is.null(newnames[[1]]) && length(newnames[[1]]) != nrow(values) || 
			!is.null(newnames[[2]]) && length(newnames[[2]]) != ncol(values)) {
			entity <- getEntityType(object)
			msg <- paste("The", entity, omxQuotes(object@name), 
				"has specified dimnames with dimensions",
				length(newnames[[1]]), "x", length(newnames[[2]]), "but the", entity,
					"is of dimensions", nrow(values), "x", ncol(values))
			stop(msg, call. = FALSE)
		}
	}
	dimnames(values) <- dimnames(object)
	return(values)
}

imxGenerateLabels <- function(model) {
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
			if(!imxIsDefinitionVariable(select[[i]]) && !hasSquareBrackets(select[[i]])) {
				labelsData[select[[i]], "model"] <- modelName
				labelsData[select[[i]], "matrix"] <- matrix@name
				labelsData[select[[i]], "row"] <- rows[[i]]
				labelsData[select[[i]], "col"] <- cols[[i]]
			}
		}
	}
	return(labelsData)
}
