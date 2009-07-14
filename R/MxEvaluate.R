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

mxEvaluate <- function(expression, model) {
	formula <- match.call()$expression
	formula <- evaluateTranslation(formula, model)
	return(eval(formula))
}

evaluateTranslation <- function(formula, model) {
	len <- length(formula)
	if (len == 0) {
		error("mxEvaluate has reached an invalid state")
	} else if (len == 1) {
		formula <- translateSymbol(formula, model)
	} else {
		formula[-1] <- lapply(formula[-1], 
			evaluateTranslation, model)
	}
	return(formula)
}

translateSymbol <- function(symbol, model) {
	key <- deparse(symbol)
	lookup <- model[[key]]
	if (is.null(lookup)) {
		return(symbol)
	} else if (is(lookup, "MxMatrix")) {
		return(lookup@values)
	} else if (is(lookup, "MxAlgebra")) {
		return(lookup@result)
	} else {
		error(paste("Cannot evaluate the object",
			omxQuotes(key), "in the model",
			omxQuotes(model@name)))
	}
}