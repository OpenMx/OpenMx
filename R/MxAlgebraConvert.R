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

convertAlgebras <- function(flatModel, convertArguments) {
    algebras <- flatModel@algebras
    if (length(algebras) == 0) {
        return(flatModel)
    }
    for(i in 1:length(algebras)) {
        flatModel <- convertSingleAlgebra(algebras[[i]], flatModel, convertArguments)
    }
    flatModel@matrices <- c(flatModel@matrices, flatModel@constMatrices)
    names(flatModel@freeMatrices) <- lapply(flatModel@freeMatrices, function(x) { x@name })
    names(flatModel@outsideMatrices) <- lapply(flatModel@outsideMatrices, function(x) { x@name })    
    flatModel@matrices <- c(flatModel@matrices, flatModel@freeMatrices, flatModel@outsideMatrices)
    return(flatModel)
}

convertSingleAlgebra <- function(algebra, flatModel, convertArguments) {
    flatModel <- convertFormulaInsertModel(algebra@formula, flatModel, convertArguments)
    formula <- convertFormula(algebra@formula, flatModel, convertArguments)
    flatModel[[algebra@name]]@formula <- formula
    return(flatModel)   
}

convertFormulaInsertModel <- function(formula, flatModel, convertArguments) {
	if (length(formula) == 1) {
        charFormula <- as.character(formula)
        if (is.numeric(formula)) {
            flatModel <- insertNumericValue(formula, flatModel)
        } else if (identical(charFormula, "")) {
            flatModel <- insertNumericValue(matrix(0,0,0), flatModel)
        } else if (charFormula %in% convertArguments$values) {
            flatModel <- insertFixedValue(charFormula, convertArguments$startvals, flatModel)
        } else if (charFormula %in% convertArguments$parameters) {
            flatModel <- insertFreeParameter(charFormula, convertArguments$startvals, flatModel)
		} else if (omxIsDefinitionVariable(charFormula)) {
             flatModel <- insertDefinitionVariable(charFormula, flatModel)
        } else if (!is.null(flatModel[[charFormula]])) {
             # do not translate this symbol
        } else if (exists(charFormula, envir = globalenv()) && is.numeric(get(charFormula, envir = globalenv()))) {
            flatModel <- insertOutsideValue(charFormula, flatModel)
        }
	} else if (length(formula) == 4 && identical(as.character(formula[1]), '[')) {
		formula[[3]] <- translateSquareBracketArgument(formula[[3]], formula[[2]], flatModel, 1)
		formula[[4]] <- translateSquareBracketArgument(formula[[4]], formula[[2]], flatModel, 2)
		for (i in 2:length(formula)) {
			flatModel <- convertFormulaInsertModel(formula[[i]], flatModel, convertArguments)
		}
	} else {
		for (i in 2:length(formula)) {
			flatModel <- convertFormulaInsertModel(formula[[i]], flatModel, convertArguments)
		}
	}
    return(flatModel)
}

insertFixedValue <- function(valName, startvals, flatModel) {
    value <- startvals[[valName]]
    flatModel <- insertNumericValue(value, flatModel)
    return(flatModel)
}

insertFreeParameter <- function(paramName, startvals, flatModel) {
    value <- as.matrix(startvals[[paramName]])
    if (!(paramName %in% names(flatModel@freeMatrices))) {
        localName <- omxUntitledName()
        identifier <- omxIdentifier(flatModel@name, localName)
        matrix <- mxMatrix("Full", values = value, labels = paramName,
            free = TRUE, name = localName)
        matrix@name <- identifier
        flatModel@freeMatrices[[paramName]] <- matrix
    }
    return(flatModel)
}

insertDefinitionVariable <- function(defName, flatModel) {
    value <- as.matrix(0)
    if (!(defName %in% names(flatModel@freeMatrices))) {
        localName <- omxUntitledName()
        identifier <- omxIdentifier(flatModel@name, localName)
        matrix <- mxMatrix("Full", values = value, labels = defName,
            free = FALSE, name = localName)
        matrix@name <- identifier
        flatModel@freeMatrices[[defName]] <- matrix
    }
    return(flatModel)
}

translateSquareBracketArgument <- function(arg, matrixName, model, rowCol) {
	if (is.character(arg)) {
		arg <- translateRowColName(matrixName, arg, model, rowCol)
	} else if (is.symbol(arg)) {
		charFormula <- as.character(arg)
		if (!identical(charFormula, '') && exists(charFormula, envir = globalenv())) {
			target <- get(charFormula, envir = globalenv())
			if (is.character(target)) {
				arg <- translateRowColName(matrixName, target, model, rowCol)
			}
		}
	}
	return(arg)
}


convertFormula <- function(formula, flatModel, convertArguments) {
	if (length(formula) == 1) {
        formula <- lookupNumericValue(formula, flatModel, convertArguments)
	} else if (length(formula) == 4 && identical(as.character(formula[1]), '[')) {
		formula[[3]] <- translateSquareBracketArgument(formula[[3]], formula[[2]], flatModel, 1)
		formula[[4]] <- translateSquareBracketArgument(formula[[4]], formula[[2]], flatModel, 2)
		for (i in 2:length(formula)) {
			formula[[i]] <- convertFormula(formula[[i]], flatModel, convertArguments)
		}	
	} else {
		for (i in 2:length(formula)) {
			formula[[i]] <- convertFormula(formula[[i]], flatModel, convertArguments)
		}
	}
    return(formula)
}

translateRowColName <- function(symbol, argname, model, rowcol) {
	key <- deparse(symbol)
	lookup <- model[[key]]
	lookupNames <- dimnames(lookup)
	if (rowcol == 1) displayRowcol <- 'row'
	else if (rowcol == 2) displayRowcol <- 'column'
	if (is.null(lookupNames)) {
		msg <- paste("The matrix", omxQuotes(key), 
			"does not contain dimnames and",
			"some algebra is referring",
			"to", displayRowcol, omxQuotes(rowcol), "by name")
		stop(msg, call. = FALSE)
	}
	if (is.null(lookupNames[[rowcol]])) {
		msg <- paste("The matrix", omxQuotes(key), 
			"does not contain", displayRowcol, "names and",
			"some algebra is referring",
			"to", displayRowcol, omxQuotes(rowcol), "by name")
		stop(msg, call. = FALSE)
	}
	rcNames <- lookupNames[[rowcol]]
	index <- match(argname, rcNames)[[1]]
	if (is.na(index)) {
		msg <- paste("The matrix", omxQuotes(key),
			"does not contain the", displayRowcol, "name",
			omxQuotes(argname))
		stop(msg, call. = FALSE)
	}
	return(index)
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

insertOutsideValue <- function(varname, flatModel) {
    value <- as.matrix(get(varname, envir = globalenv()))
    localName <- omxUntitledName()
    identifier <- omxIdentifier(flatModel@name, localName)
    matrix <- mxMatrix("Full", values = value, name = localName)
    matrix@name <- identifier
    flatModel@constMatrices[[varname]] <- matrix
    return(flatModel)
}

lookupNumericValue <- function(value, flatModel, convertArguments) {
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
    } else if (identical(as.character(value), "")) {
        for (i in 1:length(flatModel@constMatrices)) {
            constMatrix <- flatModel@constMatrices[[i]]@values
            if (nrow(constMatrix) == 0 && ncol(constMatrix) == 0) {
                return(as.symbol(flatModel@constMatrices[[i]]@name))
            }
        }
    } else if (as.character(value) %in% convertArguments$values) {
        value <- as.matrix(convertArguments$startvals[[as.character(value)]])
        for (i in 1:length(flatModel@constMatrices)) {
            constMatrix <- flatModel@constMatrices[[i]]@values
            if (nrow(value) == nrow(constMatrix) &&
                ncol(value) == ncol(constMatrix) &&
                all(value == constMatrix)) {
                return(as.symbol(flatModel@constMatrices[[i]]@name))
            }
        }
    } else if (as.character(value) %in% convertArguments$parameters) {
        matrix <- flatModel@freeMatrices[[as.character(value)]]
        return(as.symbol(matrix@name))
	} else if (omxIsDefinitionVariable(as.character(value))) {
        matrix <- flatModel@freeMatrices[[as.character(value)]]
        return(as.symbol(matrix@name))
    } else if (!is.null(flatModel[[as.character(value)]])) {
    	return(value)
    } else if (as.character(value) %in% names(flatModel@outsideMatrices)) {
        matrix <- flatModel@outsideMatrices[[as.character(value)]]
        return(as.symbol(matrix@name))    
    } else {
        return(value)
    }
}
