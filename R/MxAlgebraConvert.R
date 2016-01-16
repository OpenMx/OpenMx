#
#   Copyright 2007-2016 The OpenMx Project
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

eliminateObjectiveFunctions <- function(model) {
	model@algebras <- lapply(model@algebras, algebraEliminateObjectiveFunctions)
	if (length(model@submodels) > 0) {
		model@submodels <- lapply(model@submodels, eliminateObjectiveFunctions)
	}
	return(model)
}

algebraEliminateObjectiveFunctions <- function(algebra) {
	algebra@formula <- formulaEliminateObjectiveFunctions(algebra@formula)
	return(algebra)
}

formulaEliminateObjectiveFunctions <- function(formula) {
	if (length(formula) == 1) {
		asCharacter <- as.character(formula)
		components <- unlist(strsplit(asCharacter, imxSeparatorChar, fixed = TRUE))
		if (length(components) == 1) {
			if (identical(asCharacter, "objective")) {
				return(as.symbol("fitfunction"))
			}
		} else if (length(components) == 2) {
			modelname <- components[[1]]
			localname <- components[[2]]
			if (identical(localname, "objective")) {
				return(as.symbol(paste(modelname, "fitfunction", sep = ".")))
			}
		}
	} else {
		for (i in 2:length(formula)) {
			formula[[i]] <- formulaEliminateObjectiveFunctions(formula[[i]])
		}
	}
	return(formula)
}

constraintsToAlgebras <- function(flatModel) {
	constraints <- flatModel@constraints
	if (length(constraints) == 0) {
		return(flatModel)
	}
	for(i in 1:length(constraints)) {
		flatModel <- constraintsToAlgebrasHelper(constraints[[i]], flatModel)
	}
	return(flatModel)
}

convertAlgebras <- function(flatModel, convertArguments) {
    algebras <- flatModel@algebras
    if (length(algebras) == 0) {
        return(flatModel)
    }
    for(i in 1:length(algebras)) {
        flatModel <- convertSingleAlgebra(algebras[[i]], flatModel, convertArguments)
    }
    names(flatModel@constMatrices) <- lapply(flatModel@constMatrices, slot, "name")
    names(flatModel@freeMatrices) <- lapply(flatModel@freeMatrices, slot, "name")
    flatModel@matrices <- c(flatModel@matrices, flatModel@constMatrices, flatModel@freeMatrices)
    return(flatModel)
}

constraintsToAlgebrasHelper <- function(constraint, flatModel) {
	leftHandSide <- constraint@formula[[2]]
	rightHandSide <- constraint@formula[[3]]
	leftHandName <- imxUntitledName()
	rightHandName <- imxUntitledName()
	leftHandAlgebra <- new("MxAlgebra", NA, leftHandName, FALSE, NA_character_, NA_character_)
	rightHandAlgebra <- new("MxAlgebra", NA, rightHandName, FALSE, NA_character_, NA_character_)
	leftHandAlgebra@formula <- leftHandSide
	rightHandAlgebra@formula <- rightHandSide
	constraint@alg1 <- paste(flatModel@name, leftHandName, sep='.')
	constraint@alg2 <- paste(flatModel@name, rightHandName, sep='.')
	constraint@relation <- as.character(constraint@formula[[1]])
	flatModel[[constraint@name]] <- constraint
	flatModel[[leftHandName]] <- leftHandAlgebra
	flatModel[[rightHandName]] <- rightHandAlgebra
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
		if (!is.null(flatModel[[charFormula]])) {
             # do not translate this symbol
        } else if (is.numeric(formula)) {
            flatModel <- insertNumericValue(formula, flatModel)
        } else if (identical(charFormula, "")) {
            flatModel <- insertNumericValue(matrix(0,0,0), flatModel)
        } else if (charFormula %in% convertArguments$values) {
            flatModel <- insertFixedValue(charFormula, convertArguments$startvals, flatModel)
        } else if (charFormula %in% convertArguments$parameters) {
            flatModel <- insertFreeParameter(charFormula, convertArguments$startvals,
		convertArguments$bounds, flatModel)
		} else if (imxIsDefinitionVariable(charFormula)) {
             flatModel <- insertDefinitionVariable(charFormula, flatModel)
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

insertFreeParameter <- function(paramName, startvals, bounds, flatModel) {
    if (!(paramName %in% names(flatModel@freeMatrices))) {
		localName <- imxUntitledName()
		identifier <- imxIdentifier(flatModel@name, localName)
		value <- as.matrix(startvals[[paramName]])
		lbound <- as.matrix(bounds[[paramName]][[1]])
		ubound <- as.matrix(bounds[[paramName]][[2]])
        matrix <- mxMatrix("Full", values = value, labels = paramName,
            free = TRUE, lbound = lbound, ubound = ubound, name = localName)
        matrix@name <- identifier
		matrix@display <- paramName
		matrix <- findSquareBrackets(matrix)
		flatModel@freeMatrices[[paramName]] <- matrix
    }
    return(flatModel)
}

insertDefinitionVariable <- function(defName, flatModel) {
    if (!(defName %in% names(flatModel@freeMatrices))) {
        localName <- imxUntitledName()
        identifier <- imxIdentifier(flatModel@name, localName)
	    value <- as.matrix(0)
        matrix <- mxMatrix("Full", values = value, labels = defName,
            free = FALSE, name = localName)
        matrix@name <- identifier
		matrix@display <- defName
		matrix <- findSquareBrackets(matrix)
        flatModel@freeMatrices[[defName]] <- matrix
    }
    return(flatModel)
}

squareBracketArgumentHelper <- function(arg, matrixName, model, rowCol) {
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

translateSquareBracketArgument <- function(arg, matrixName, model, rowCol) {
	if (length(arg) > 1) {
		for (i in 1:length(arg)) {
			arg[[i]] <- translateSquareBracketArgument(arg[[i]], matrixName, model, rowCol)
		}
		return(arg)
	} else {
		return(squareBracketArgumentHelper(arg, matrixName, model, rowCol))
	}
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
    if (length(flatModel@constMatrices) == 0) {
        localName <- imxUntitledName()
        identifier <- imxIdentifier(flatModel@name, localName)
        matrix <- mxMatrix("Full", values = as.matrix(value), name = localName)
        matrix@name <- identifier
		matrix@display <- as.character(value)
		matrix <- findSquareBrackets(matrix)
        flatModel@constMatrices[[identifier]] <- matrix
    } else {
		valuematrix <- as.matrix(value)
        for (i in 1:length(flatModel@constMatrices)) {
            constMatrix <- flatModel@constMatrices[[i]]@values
            if (nrow(valuematrix) == nrow(constMatrix) &&
                ncol(valuematrix) == ncol(constMatrix) &&
                all(valuematrix == constMatrix)) {
                return(flatModel)
            }
        }
        localName <- imxUntitledName()
        identifier <- imxIdentifier(flatModel@name, localName)
        matrix <- mxMatrix("Full", values = valuematrix, name = localName)
        matrix@name <- identifier
		matrix@display <- as.character(value)
		matrix <- findSquareBrackets(matrix)
        flatModel@constMatrices[[identifier]] <- matrix
    }
    return(flatModel)
}

insertOutsideValue <- function(varname, flatModel) {
    value <- as.matrix(get(varname, envir = globalenv()))
    if (length(flatModel@constMatrices) == 0) {
	    localName <- imxUntitledName()
	    identifier <- imxIdentifier(flatModel@name, localName)
        matrix <- mxMatrix("Full", values = value, name = localName)
	    matrix@name <- identifier
		matrix@display <- varname
		matrix <- findSquareBrackets(matrix)
    	flatModel@constMatrices[[varname]] <- matrix
    } else {
        for (i in 1:length(flatModel@constMatrices)) {
            constMatrix <- flatModel@constMatrices[[i]]@values
            if (nrow(value) == nrow(constMatrix) &&
                ncol(value) == ncol(constMatrix) &&
                all(value == constMatrix)) {
                return(flatModel)
            }
        }
        localName <- imxUntitledName()
	    identifier <- imxIdentifier(flatModel@name, localName)
        matrix <- mxMatrix("Full", values = value, name = localName)
	    matrix@name <- identifier
		matrix@display <- varname
		matrix <- findSquareBrackets(matrix)
    	flatModel@constMatrices[[varname]] <- matrix
    }
    return(flatModel)
}

lookupNumericValue <- function(value, flatModel, convertArguments) {
   asCharacter <- as.character(value)
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
    } else if (identical(asCharacter, "")) {
        for (i in 1:length(flatModel@constMatrices)) {
            constMatrix <- flatModel@constMatrices[[i]]@values
            if (nrow(constMatrix) == 0 && ncol(constMatrix) == 0) {
                return(as.symbol(flatModel@constMatrices[[i]]@name))
            }
        }
    } else if (asCharacter %in% convertArguments$values) {
        value <- as.matrix(convertArguments$startvals[[asCharacter]])
        for (i in 1:length(flatModel@constMatrices)) {
            constMatrix <- flatModel@constMatrices[[i]]@values
            if (nrow(value) == nrow(constMatrix) &&
                ncol(value) == ncol(constMatrix) &&
                all(value == constMatrix)) {
                return(as.symbol(flatModel@constMatrices[[i]]@name))
            }
        }
    } else if (asCharacter %in% convertArguments$parameters) {
        matrix <- flatModel@freeMatrices[[asCharacter]]
        return(as.symbol(matrix@name))
	} else if (imxIsDefinitionVariable(asCharacter)) {
        matrix <- flatModel@freeMatrices[[asCharacter]]
        return(as.symbol(matrix@name))
	} else if(exists(asCharacter, envir = globalenv()) && 
		is.numeric(get(asCharacter, envir = globalenv()))) {
		value <- as.matrix(get(asCharacter, envir = globalenv()))
        for (i in 1:length(flatModel@constMatrices)) {
            constMatrix <- flatModel@constMatrices[[i]]@values
            if (nrow(value) == nrow(constMatrix) &&
                ncol(value) == ncol(constMatrix) &&
                all(value == constMatrix)) {
                return(as.symbol(flatModel@constMatrices[[i]]@name))
            }
        }
    }
    return(value)
}
