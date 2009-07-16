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
	formula <- evaluateTranslation(formula, model, modelVariable)
	if (show) {
		cat(deparse(formula, width.cutoff = 500L), '\n')
	}
	return(eval(formula))
}

evaluateTranslation <- function(formula, model, modelVariable) {
	len <- length(formula)
	if (len == 0) {
		stop("mxEvaluate has reached an invalid state")
	} else if (len == 1) {
		formula <- translateSymbol(formula, model, modelVariable)
	} else {
		formula[-1] <- lapply(formula[-1], 
			evaluateTranslation, model, modelVariable)
	}
	return(formula)
}

translateSymbol <- function(symbol, model, modelVariable) {
	key <- deparse(symbol)
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
