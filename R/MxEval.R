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

mxEval <- function(expression, model, compute = FALSE, show = FALSE, defvar.row = 1L,
			cache = new.env(parent = emptyenv()), cacheBack = FALSE, .extraBack=0L) {
	if (missing(expression)) {
		stop("'expression' argument is mandatory in call to mxEval function")
	} else if (missing(model)) {
		stop("'model' argument is mandatory in call to mxEval function")
	}
	expression <- match.call()$expression
	modelvariable <- match.call()$model
  #for (x in 1:sys.nframe()) cat(x-1, head(ls(parent.frame(x))), fill=TRUE)
	EvalInternal(expression, model, modelvariable, compute, show, defvar.row,
		       cache, cacheBack, 2L + .extraBack)
}

EvalInternal <- function(expression, model, modelvariable, compute, show, defvar.row,
		       cache, cacheBack, back) {
	if (!is.numeric(defvar.row) || length(defvar.row) != 1 || is.na(defvar.row)) {
		stop("'defvar.row' argument must be a single numeric value")
	}
	labelsData <- imxGenerateLabels(model)
	env <- parent.frame(back)
	if (compute) {
		model <- imxPreprocessModel(model)
		expression <- formulaEliminateObjectiveFunctions(expression)
		namespace <- imxGenerateNamespace(model)
		model <- imxFlattenModel(model, namespace, FALSE)
		expression <- qualifyNamesFormula(expression, model@name, namespace)
	}
	if (show) {
		tuple <- evaluateExpression(expression, deparse(expression), model, 
			labelsData, env, compute, show = TRUE, outsideAlgebra = TRUE,
			cache = new.env(parent = emptyenv()))
		showResult <- tuple[[1]]
		showResult <- eval(substitute(substitute(x, list(.zzz = modelvariable)), list(x = showResult)))
		print(deparse(showResult), width.cutoff = 450L)
	}
	tuple <- evaluateExpression(expression, deparse(expression), model, 
		labelsData, env, compute, show = FALSE, outsideAlgebra = TRUE, defvar.row,
		cache = cache)
	result <- tuple[[1]]
	cache <- tuple[[2]]
	if (!is.vector(result)) {
		result <- eval(result, envir = env)
	}
	if (!is.matrix(result)) {
    if (is(result, "MxAlgebra")) {
      stop(paste("You have referenced algebra",
                 omxQuotes(result$name), "outside of the model"))
    }
    result <- as.matrix(result)
  }
	if (cacheBack) {
		return(list(result, cache))
	} else {
		return(result)
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
		if(!identical(as.character(formula), "")){
			formula <- translateErrorSymbol(formula, model)
		}
	} else {
		formula <- lapply(formula, translateErrorFormula, model)
		formula <- as.call(formula)
	}
	return(formula)
}

evaluateExpression <- function(formula, contextString, model, labelsData, 
	env, compute, show, outsideAlgebra, defvar.row = 1, cache = new.env(parent = emptyenv())) {
	len <- length(formula)
	if (len == 0) {
		stop("mxEval has reached an invalid state")
	} else if (len == 1) {
		if (!identical(as.character(formula), "")) {
			tuple <- evaluateSymbol(formula, contextString, model, 
				labelsData, env, compute, show, outsideAlgebra, defvar.row,
				cache)
			return(tuple)
		} else {
			return(list(formula, cache))
		}
	}
	originalFormula <- formula
	for(i in 2:length(formula)) {
		tuple <- evaluateExpression(formula[[i]], contextString, 
			model, labelsData, env, compute, show, 
			outsideAlgebra, defvar.row, cache)
		formula[[i]] <- tuple[[1]]
		cache <- tuple[[2]]
	}
	if (len == 4 && identical(as.character(formula[1]), '[')) {
		formula$drop <- FALSE
		if(is.matrix(formula[[3]]) && nrow(formula[[3]]) == 0 && ncol(formula[[3]]) == 0) {
			formula[[3]] <- TRUE
		}
		if(is.matrix(formula[[4]]) && nrow(formula[[4]]) == 0 && ncol(formula[[4]]) == 0) {
			formula[[4]] <- TRUE
		}
	}
	operator <- as.character(formula[1])
	if (len == 3 && operator %in% c('*', '^', '/', '+', '-')) {
		if (!is.numeric(formula[[2]]) || length(formula[[2]]) > 1){
			formula[[2]] <- substitute(as.matrix(x), list(x = formula[[2]]))
		} else if (is.matrix(formula[[2]]) && length(formula[[2]]) == 1) {
			formula[[2]] <- substitute(as.numeric(x), list(x = formula[[2]]))
		}
		if (!is.numeric(formula[[3]]) || length(formula[[3]]) > 1){
			formula[[3]] <- substitute(as.matrix(x), list(x = formula[[3]]))
		} else if (is.matrix(formula[[3]]) && length(formula[[3]]) == 1) {
			formula[[3]] <- substitute(as.numeric(x), list(x = formula[[3]]))
		}
	}
	if (show) {
		return(list(result, cache))
	}
	result <- tryCatch(do.call(as.character(formula[[1]]), as.list(formula[-1]), envir = env),
			error = function(x) {
				originalFormula <- translateErrorFormula(originalFormula, model)
				formulaString <- deparse(originalFormula, width.cutoff = 500L)
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
	return(list(result, cache))
}

evaluateSymbol <- function(symbol, contextString, model, labelsData,
			env, compute, show, outsideAlgebra, defvar.row = 1, 
			cache) {
	key <- deparse(symbol)
	if (exists(key, envir = cache)) {
		result <- get(key, envir = cache)
	} else {
		index <- match(key, dimnames(labelsData)[[1]])
		if (!is.na(index)) {
			targetmodel <- labelsData[index, "model"]
			matrix <- labelsData[index, "matrix"]
			fullname <- imxIdentifier(targetmodel, matrix)
			row <- labelsData[index, "row"]
			col <- labelsData[index, "col"]
			value <- model[[fullname]]@values[row,col]
			if (show) {
				result <- substitute(.zzz[[x]]@values[y,z], list(x = fullname, y = row, z = col))
			} else {
				result <- value
			}
		} else {
			lookup <- model[[key]]
			if (is.list(lookup) && length(lookup) == 2 && is(lookup[[2]], "MxFitFunction")) {
				lookup <- lookup[[2]]
			}
			if (is.list(lookup) && length(lookup) == 2 && is(lookup[[2]], "MxExpectation")) {
				lookup <- lookup[[2]]
			}
			if (imxIsDefinitionVariable(key)) {
				result <- definitionStartingValue(key, contextString, model, defvar.row)
			} else if (is.null(lookup)) {
				if (!show && !outsideAlgebra && exists(key, envir = env)) {
					result <- try(as.matrix(get(key, envir = env)), silent=TRUE)
          if (is(result, 'try-error')){
          	stop(paste(
          		"cannot coerce object with symbol", omxQuotes(key), "to a matrix\n",
          		"(hint: did you neglect to put an object named", omxQuotes(key), "into model",
          		omxQuotes(model@name),"?)"))
          }
				} else {
					result <- symbol
				}
			} else if (is(lookup, "MxMatrix")) {
				tuple <- evaluateMatrix(lookup, model, labelsData, env, show, compute, defvar.row, cache)
				result <- tuple[[1]]
				cache <- tuple[[2]]
			} else if (is(lookup, "MxAlgebra")) {
				tuple <- evaluateAlgebra(lookup, model, labelsData, env, show, compute, defvar.row, cache)
				result <- tuple[[1]]
				cache <- tuple[[2]]
			} else if (is(lookup, "MxData")) {
				if (show) {
					result <- substitute(.zzz[[x]]@observed, list(x = key))
				} else {
					result <- lookup@observed
				}
			} else if (is(lookup, "MxFitFunction")) {
				if (length(lookup@result) == 0) {
					if (compute) {
						result <- genericFitInitialMatrix(lookup, model)
					} else {
						result <- matrix(0, 0, 0)
					}
				} else if (show) {
					result <- substitute(.zzz[[x]]@result, list(x = key))
				} else {
					result <- lookup@result
				}
			} else if (is(lookup, "MxExpectation")) {
				stop("mxEval for Expectations is not yet implemented.")
				#signpost
				#if (length(lookup@result) == 0) {
				#	if (compute) {
				#		result <- genericFitInitialMatrix(lookup, model)
				#	} else {
				#		result <- matrix(0, 0, 0)
				#	}
				#} else if (show) {
				#	result <- substitute(.zzz[[x]]@result, list(x = key))
				#} else {
				#	result <- lookup@result
				#}
			} else {
				stop(paste("Cannot evaluate the object",
					omxQuotes(key), "in the model",
					omxQuotes(model@name)))
			}
		}
		assign(key, result, envir = cache)
	}
	return(list(result, cache))
}

evaluateMatrix <- function(matrix, model, labelsData, env, show, compute, defvar.row, cache) {
	if (show) {
		return(list(substitute(.zzz[[x]]@values, list(x = matrix@name)), cache))
	} else if (compute) {
		tuple <- computeMatrix(matrix, model, labelsData, defvar.row, env, cache)
	} else {
		tuple <- list(matrix@values, cache)
	}
	tuple[[1]] <- assignDimnames(matrix, tuple[[1]])
	return(tuple)
}

computeMatrix <- function(matrix, model, labelsData, defvar.row, env, cache) {
	matrix <- populateDefVarMatrix(matrix, model, defvar.row)
	values <- matrix@values
	labels <- matrix@labels
	select <- matrix@.squareBrackets
	if (all(!select)) {
		return(list(values, cache))
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
		tuple <- evaluateExpression(expression, contextString, model, labelsData, env, 
			compute=TRUE, show=FALSE, outsideAlgebra=FALSE, defvar.row=defvar.row, cache=cache)
		result <- tuple[[1]]
		cache <- tuple[[2]]
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
	return(list(values, cache))
}

evaluateAlgebra <- function(algebra, model, labelsData, env, show, compute, defvar.row, cache) {
	if (compute && show) {
		return(computeAlgebra(algebra, model, labelsData, show = TRUE, defvar.row, env, cache))
	} else if (compute) {
		tuple <- computeAlgebra(algebra, model, labelsData, show = FALSE, defvar.row, env, cache)
		tuple[[1]] <- as.matrix(tuple[[1]])
	} else if (show) {
		result <- substitute(.zzz[[x]]@result, list(x = algebra@name))
		return(list(result, cache))
	} else {
		tuple <- list(algebra@result, cache)
	}
	if (!is.null(dimnames(algebra))) {
		tuple[[1]] <- assignDimnames(algebra, tuple[[1]])
	}
	return(tuple)
}

computeAlgebra <- function(algebra, model, labelsData, show, defvar.row, env, cache) {
	contextString <- simplifyName(algebra@name, model@name)
	tuple <- evaluateExpression(algebra@formula, contextString, model, labelsData, env, 
		compute = TRUE, show, outsideAlgebra = FALSE, defvar.row, cache)
	return(tuple)
}

evaluateMxObject <- function(objname, flatModel, labelsData, cache) {
	if (flatModel@unsafe) {
		#stop("evaluateMxObject must be avoided in unsafe context")
		# mxFitFunctionRow genericFitAddEntities still needs this
	}
	return(eval(substitute(evaluateSymbol(x, objname, flatModel, 
			labelsData, globalenv(), compute = TRUE, 
			show = FALSE, outsideAlgebra = FALSE, defvar.row = 1, 
			cache = cache),
			list(x = quote(as.symbol(objname))))))
}

evaluateAlgebraWithContext <- function(algebra, context, flatModel, labelsData, cache) {
	return(evaluateExpression(algebra@formula, context, flatModel, labelsData, globalenv(), 
		compute = TRUE, show = FALSE, outsideAlgebra = FALSE, cache = cache))
}

assignDimnames <- function(object, values) {
	newnames <- dimnames(object)
	if (!is.null(newnames)) {
		if (length(newnames) != 2) {
			entity <- class(object)
			msg <- paste("The 'dimnames' argument to", entity,
				omxQuotes(object@name), "must have a length of 2")
			stop(msg, call. = FALSE)
		}
		if (!is.null(newnames[[1]]) && length(newnames[[1]]) != nrow(values) || 
			!is.null(newnames[[2]]) && length(newnames[[2]]) != ncol(values)) {
			entity <- class(object)
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

##' imxGenerateLabels
##'
##' This is an internal function exported for those people who know
##' what they are doing.
##'
##' @param model model
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

retainLabel <- function(label, existing) {
    return(!imxIsDefinitionVariable(label) && 
			!hasSquareBrackets(label) &&
			!(label %in% existing))
}

generateLabelsMatrix <- function(modelName, matrix, labelsData) {
	labels <- matrix@labels
	test <- !is.na(labels)
	select <- labels[test]
	rows <- row(labels)[test]
	cols <- col(labels)[test]
	if (length(select) > 0) {
		retain <- !duplicated(select) & sapply(select, 
			retainLabel, rownames(labelsData))
		select <- select[retain]
		rows <- rows[retain]
		cols <- cols[retain]
		if (length(select) > 0) {
			newLabels <- data.frame(model=modelName, 
				matrix=matrix@name, row=rows, 
				col=cols, stringsAsFactors=FALSE)
			rownames(newLabels) <- select
			labelsData <- rbind(labelsData, newLabels)
		}
	}
	return(labelsData)
}

mxEvalByName <- function(name, model, compute=FALSE, show=FALSE, defvar.row = 1L,
		cache = new.env(parent = emptyenv()), cacheBack = FALSE, .extraBack=0L) {
   if((length(name) != 1) || typeof(name) != "character") {
      stop("'name' argument must be a character argument")
   }
   warnModelCreatedByOldVersion(model)
#for (x in 1:sys.nframe()) cat(x-1, head(ls(parent.frame(x))), fill=TRUE)
   eval(substitute(mxEval(x, model, compute, show, defvar.row, cache, cacheBack, 1L + .extraBack),
      list(x = as.symbol(name))))
}

BootstrapEvalInternal <- function(expression, model, modelvariable, defvar.row, ...) {
	bootData <- omxGetBootstrapReplications(model)
  mle <- cvectorize(EvalInternal(expression, model, modelvariable,
				 compute=TRUE, show=FALSE, defvar.row, cache=new.env(parent = emptyenv()),
						 cacheBack=FALSE, back=2L))
  result <- matrix(NA, nrow=nrow(bootData), ncol=length(mle))
  for (rx in 1:nrow(bootData)) {
	  bmodel <- omxSetParameters(model, labels=colnames(bootData), values=unlist(bootData[rx,]))
	  result[rx,] <- cvectorize(EvalInternal(expression, bmodel, modelvariable,
						 compute=TRUE, show=FALSE, defvar.row, cache=new.env(parent = emptyenv()),
						 cacheBack=FALSE, back=2L))
  }
	list(mle=mle, result=result)
}

omxBootstrapEval <- function(expression, model, defvar.row = 1L, ...) {
  expression <- match.call()$expression
  modelvariable <- match.call()$model
  BootstrapEvalInternal(expression, model, modelvariable, defvar.row, ...)$result
}

omxBootstrapEvalCov <- function(expression, model, defvar.row = 1L, ...) {
  expression <- match.call()$expression
  modelvariable <- match.call()$model
  cov(BootstrapEvalInternal(expression, model, modelvariable, defvar.row, ...)$result)
}

mxBootstrapEval <- function(expression, model, defvar.row = 1L, ..., bq=c(.25,.75),
				method=c('bcbci','quantile')) {
  warnModelCreatedByOldVersion(model)
  expression <- match.call()$expression
  modelvariable <- match.call()$model
	method <- match.arg(method)
  got <- BootstrapEvalInternal(expression, model, modelvariable, defvar.row, ...)
	cbind("SE"=apply(got$result, 2, sd),
	      summarizeBootstrap(got$mle, got$result, bq, method))
}

mxBootstrapEvalByName <- function(name, model, defvar.row=1L, ..., bq=c(.25,.75),
				  method=c('bcbci','quantile')) {
   if((length(name) != 1) || typeof(name) != "character") {
      stop("'name' argument must be a character argument")
   }
  warnModelCreatedByOldVersion(model)
   method <- match.arg(method)
   eval(substitute(mxBootstrapEval(x, model, defvar.row, bq=bq, method=method),
      list(x = as.symbol(name))))
}

omxBootstrapEvalByName <- function(name, model, defvar.row=1L, ...) {
   if((length(name) != 1) || typeof(name) != "character") {
      stop("'name' argument must be a character argument")
   }
   method <- match.arg(method)
   eval(substitute(omxBootstrapEval(x, model, defvar.row),
      list(x = as.symbol(name))))
}
