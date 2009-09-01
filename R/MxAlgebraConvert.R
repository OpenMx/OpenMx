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

convertAlgebras <- function(flatModel) {
    algebras <- flatModel@algebras
    if (length(algebras) == 0) {
        return(flatModel)
    }
    for(i in 1:length(algebras)) {
        flatModel <- convertSingleAlgebra(algebras[[i]], flatModel)
    }
    flatModel@matrices <- c(flatModel@matrices, flatModel@constMatrices)
    return(flatModel)
}

convertSingleAlgebra <- function(algebra, flatModel) {
    flatModel <- convertFormulaInsertModel(algebra@formula, flatModel)
    formula <- convertFormula(algebra@formula, flatModel)
    flatModel[[algebra@name]]@formula <- formula
    return(flatModel)   
}

convertFormulaInsertModel <- function(formula, flatModel) {
	if (length(formula) == 1) {
        if (is.numeric(formula)) {
            flatModel <- insertNumericValue(formula, flatModel)
        }
	} else {
		for (i in 2:length(formula)) {
			flatModel <- convertFormulaInsertModel(formula[[i]], flatModel)
		}
	}
    return(flatModel)
}

convertFormula <- function(formula, flatModel) {
	if (length(formula) == 1) {
        formula <- lookupNumericValue(formula, flatModel)
	} else {
		for (i in 2:length(formula)) {
			formula[[i]] <- convertFormula(formula[[i]], flatModel)
		}
	}
    return(formula)
}

insertNumericValue <- function(value, flatModel) {
    value <- as.matrix(value)
    if (length(flatModel@constMatrices) == 0) {
        localName <- omxUntitledName()
        identifier <- omxIdentifier(flatModel@name, localName)
        matrix <- mxMatrix("Full", values = value, name = localName)
        matrix@name <- identifier
        flatModel@constMatrices[[identifier]] <- matrix
    } else {
        for (i in 1:length(flatModel@constMatrices)) {
            constMatrix <- flatModel@constMatrices[[i]]@values
            if (nrow(value) == nrow(constMatrix) &&
                ncol(value) == ncol(constMatrix) &&
                all(value == constMatrix)) {
                return(flatModel)
            }
        }
        localName <- omxUntitledName()
        identifier <- omxIdentifier(flatModel@name, localName)
        matrix <- mxMatrix("Full", values = value, name = localName)
        matrix@name <- identifier
        flatModel@constMatrices[[identifier]] <- matrix
    }
    return(flatModel)
}

lookupNumericValue <- function(value, flatModel) {
    if (is.numeric(value)) {
        value <- as.matrix(value)
        for (i in 1:length(flatModel@constMatrices)) {
            constMatrix <- flatModel@constMatrices[[i]]@values
            if (nrow(value) == nrow(constMatrix) &&
                ncol(value) == ncol(constMatrix) &&
                all(value == constMatrix)) {
                return(as.symbol(flatModel@constMatrices[[i]]@name))
            }
        }
    } else {
        return(value)
    }
}

