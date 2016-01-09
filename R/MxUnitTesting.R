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


checkZeroDimensions <- function(a, b) {
	if ((length(a) == 0 && length(b) > 0) ||
		(length(a) > 0 && length(b) == 0)) {
			stop(paste("One of these has zero length:",
				omxQuotes(paste(a, collapse = ' ')), 
				"and", omxQuotes(paste(b, collapse = ' ')))) 	
	} else if (length(a) == 0 && length(b) == 0) {
		warning("Both values have zero length.  That's weird.")
	}
}

checkEqualDimensions <- function(a, b) {
	checkZeroDimensions(a, b)
	if((is.vector(a) && length(a) > 1 && !is.vector(b)) || 
		(is.vector(b) && length(b) > 1 && !is.vector(a))) {
		stop(paste("Not both vectors:",
			   omxQuotes(paste(a, collapse = ' ')),
			   "and", omxQuotes(paste(b, collapse = ' '))))
	}
	if((is.matrix(a) && (nrow(a) > 1 || ncol(a) > 1) && !is.matrix(b)) || 
		(is.matrix(b) && (nrow(b) > 1 || ncol(b) > 1) && !is.matrix(a))) {
		stop(paste("Not both matrices:",
			   omxQuotes(paste(a, collapse = ' ')),
			   "and", omxQuotes(paste(b, collapse = ' '))))
	}	
	if (is.vector(a) && (length(a) != length(b))) {
		stop(paste("Not equal length",
			   length(a), 'and', length(b), ":",
			   omxQuotes(paste(a, collapse = ' ')),
			   "and", omxQuotes(paste(b, collapse = ' '))))
	}
	if (is.matrix(a) && (nrow(a) > 1 || ncol(a) > 1) && any(dim(a) != dim(b))) {
		stop(paste("Not of equal dimension",
			   paste(dim(a), collapse = ' x '), 'and',
			   paste(dim(b), collapse = ' x '), ":",
			   omxQuotes(paste(a, collapse = ' ')),
			   "and", omxQuotes(paste(b, collapse = ' '))))
	}
}

omxCheckIdentical <- function(a, b) {
	checkEqualDimensions(a, b)	
	if (any(!identical(a, b))) {
		stop(paste(omxQuotes(paste(a, collapse = ' ')), 
			"and", omxQuotes(paste(b, collapse = ' ')), 
			"are not identical"))
	} else if (getOption("mxPrintUnitTests")) {
		cat(paste(deparse(match.call()$a), "and", 
			deparse(match.call()$b),
			"are identical.\n"))	
	}
}

omxCheckEquals <- function(a, b) {
	checkEqualDimensions(a, b)	
	if (any(a != b)) {
		stop(paste(omxQuotes(paste(a, collapse = ' ')), 
			"and", omxQuotes(paste(b, collapse = ' ')), 
			"are not equal"))
	} else if (getOption("mxPrintUnitTests")) {
		cat(paste(deparse(match.call()$a), "and", 
			deparse(match.call()$b),
			"are equal.\n"))	
	}
}

omxCheckSetEquals <- function(a, b) {
	if (!setequal(a, b)) {
		stop(paste(omxQuotes(paste(b, collapse = ' ')), 
			"and", omxQuotes(paste(b, collapse = ' ')), 
			"do not contain the same elements"))
	} else if (getOption("mxPrintUnitTests")) {
		cat(paste(deparse(match.call()$a), "and", 
			deparse(match.call()$b),
			"contain the same elements.\n"))
	}
}

omxCheckTrue <- function(a) {	
	if (any(!a)) {
		call <- deparse(match.call()$a)
		stop(paste(call, "is not true"))
	} else if (getOption("mxPrintUnitTests")) {
		call <- deparse(match.call()$a)
		cat(paste(call, "is true.", '\n'))
	}
}

##' Approximate Equality Testing Function
##'
##' This function tests whether two numeric vectors or matrixes are
##' approximately equal to one another, within a specified threshold.
##'
##' Arguments \sQuote{a} and \sQuote{b} must be of the same type,
##' ie. they must be either vectors of equal dimension or matrices of
##' equal dimension. The two arguments are compared element-wise for
##' approximate equality.  If the absolute value of the difference of
##' any two values is greater than the threshold, then an error will
##' be thrown. If \sQuote{a} and \sQuote{b} are approximately equal to
##' each other, by default the function will print a statement
##' informing the user the test has passed.  To turn off these print
##' statements use \code{options("mxPrintUnitTests" = FALSE)}.
##'
##' When na.action is set to na.pass, a and b are expected to have
##' identical missingness patterns.
##'
##' @param a a numeric vector or matrix
##' @param b a numeric vector or matrix
##' @param epsilon a non-negative tolerance threshold
##' @param na.action either na.fail (default) or na.pass. Use of
##' na.omit or na.exclude is not recommended.
##' @seealso
##' \code{\link{omxCheckWithinPercentError}},
##' \code{\link{omxCheckIdentical}}, \code{\link{omxCheckSetEquals}},
##' \code{\link{omxCheckTrue}}, \code{\link{omxCheckEquals}}
##' @references
##' The OpenMx User's guide can be found at http://openmx.psyc.virginia.edu/documentation.
##' @examples
##' omxCheckCloseEnough(c(1, 2, 3), c(1.1, 1.9 ,3.0), epsilon = 0.5)
##' omxCheckCloseEnough(matrix(3, 3, 3), matrix(4, 3, 3), epsilon = 2)
##' # Throws an error
##' try(omxCheckCloseEnough(c(1, 2, 3), c(1.1, 1.9 ,3.0), epsilon = 0.01))
omxCheckCloseEnough <- function(a, b, epsilon = 10^(-15), na.action=na.fail) {
	if (epsilon < 0) stop("epsilon must be non-negative")
	checkEqualDimensions(a, b)
	a <- na.action(a)
	b <- na.action(b)
	if (any(is.na(a) != is.na(b))) {
		stop(paste("In", deparse(width.cutoff = 400L, sys.call()), ":",
			   "different missingness patterns:",
			   omxQuotes(paste(a, collapse = ' ')),
			   "and", omxQuotes(paste(b, collapse = ' '))),
			call. = FALSE)
	}
	a.vec <- as.vector(a)
	b.vec <- as.vector(b)
	a.vec <- a.vec[!is.na(a.vec)]
	b.vec <- b.vec[!is.na(b.vec)]
	check <- any(mapply(function(x,y) {
			abs(x - y) > epsilon }, 
			a.vec, b.vec))
	if (check) {
		stop(paste("In", deparse(width.cutoff = 400L, sys.call()), ":",
			   "not equal to within", epsilon, ':',
			   omxQuotes(paste(a, collapse = ' ')),
			   "and", omxQuotes(paste(b, collapse = ' '))),
			call. = FALSE)
	} else if (getOption("mxPrintUnitTests")) {
		cat(paste(deparse(match.call()$a), "and", 
			deparse(match.call()$b),
			"are equal to within", paste(epsilon, ".\n", sep = '')))
	}
}

omxCheckWithinPercentError <- function(a, b, percent = 0.1) {
	checkEqualDimensions(a, b)	
	check <- any(mapply(function(x,y) {
			(abs(x - y)/x * 100) > percent }, 
			as.vector(a), as.vector(b)))	
	if (check) {
		stop(paste(omxQuotes(paste(a, collapse = ' ')), 
			"does not estimate", 
			omxQuotes(paste(b, collapse = ' ')), 
			"within", percent, "percent"))
	} else if (getOption("mxPrintUnitTests")) {
		cat(paste(deparse(match.call()$a), "and", 
			deparse(match.call()$b),
			"are equal to within", percent, "percent.\n"))
	}
}

trim <- function(input) {
	input <- sub("(?m)^\\s+", "", input, perl = TRUE)
	input <- sub("(?m)\\s+$", "", input, perl = TRUE)
	return(input)
}

tryCatch.W <- function(expr) {
	# see demo(error.catching)
	W <- NULL
	w.handler <- function(w) {
		W <<- w
		invokeRestart("muffleWarning")
	}
	list(value = withCallingHandlers(tryCatch(expr), warning = w.handler),
	     warning = W)
}

##' Correct Warning Message Function
##'
##' This function tests whether the correct warning message is thrown.
##'
##' Arguments \sQuote{expression} and \sQuote{message} give the expression
##' that generates the warning and the message that is supposed to be generated, respectively.
##'
##' @param expression an R expression that produces a warning
##' @param message a character string with the desired warning message
##' @seealso
##' \code{\link{omxCheckError}}
##' \code{\link{omxCheckWithinPercentError}},
##' \code{\link{omxCheckIdentical}}, \code{\link{omxCheckSetEquals}},
##' \code{\link{omxCheckTrue}}, \code{\link{omxCheckEquals}}
##' @references
##' The OpenMx User's guide can be found at http://openmx.psyc.virginia.edu/documentation.
##' @examples
##' msg <- paste("Objective functions like mxFIMLObjective()",
##'		"have been deprecated in favor of expectation and fit functions.\n",
##'		"Please use mxExpectationNormal(covariance= , means = , ...) instead,",
##'		"and add a call to mxFitFunctionML().",
##'		"See examples at help(mxExpectationNormal)")
##' foo <- omxCheckWarning(mxFIMLObjective('cov', 'mean'), msg)
omxCheckWarning <- function(expression, message) {
	inputExpression <- match.call()$expression
	result <- tryCatch.W(expression)
	if (is.null(result$warning)) {
		if (all(is.na(message))) {
			if (getOption("mxPrintUnitTests")) {
				cat(paste("As expected, no warning was generated by the expression",
					  deparse(inputExpression, width.cutoff = 500L), '\n'))
			}
			return(result$value)
		}
		stop(paste("No warning was observed for the expression",
			deparse(inputExpression, width.cutoff = 500L)), call. = FALSE)
	}
	if (trim(result$warning$message) != trim(message)) {
		stop(paste("A warning was thrown with the wrong message:",
			   result$warning$message), call. = FALSE)
	}
	if (getOption("mxPrintUnitTests")) {
		cat(paste("The expected warning was observed for the expression",
			  deparse(inputExpression, width.cutoff = 500L), '\n'))
	}
	result$value
}

##' Correct Error Message Function
##'
##' This function tests whether the correct error message is thrown.
##'
##' Arguments \sQuote{expression} and \sQuote{message} give the expression
##' that generates the error and the message that is supposed to be generated, respectively.
##'
##' @param expression an R expression that produces an error
##' @param message a character string with the desired error message
##' @seealso
##' \code{\link{omxCheckWarning}}
##' \code{\link{omxCheckWithinPercentError}},
##' \code{\link{omxCheckIdentical}}, \code{\link{omxCheckSetEquals}},
##' \code{\link{omxCheckTrue}}, \code{\link{omxCheckEquals}}
##' @references
##' The OpenMx User's guide can be found at http://openmx.psyc.virginia.edu/documentation.
##' @examples
##' A <- mxMatrix('Full', 1, 1, labels = 'data.foo', free = TRUE, name = 'A')
##' model <- mxModel('model', A)
##' omxCheckError(mxRun(model), 
##'	paste("The definition variable 'data.foo'",
##'		"has been assigned to a",
##'		"free parameter in matrix 'A'"))
##' omxCheckCloseEnough(matrix(3, 3, 3), matrix(4, 3, 3), epsilon = 2)
##' # Throws error, check the message
##' tmsg <- paste("In omxCheckCloseEnough(c(1, 2, 3), c(1.1, 1.9, 3), 0.01)",
##'		": not equal to within 0.01 : '1 2 3' and '1.1 1.9 3'")
##' omxCheckError(omxCheckCloseEnough(c(1, 2, 3), c(1.1, 1.9 ,3.0), .01), tmsg)
omxCheckError <- function(expression, message) {
	inputExpression <- match.call()$expression
	assign("checkErrorState", FALSE, pkg_globals)
	tryCatch(eval(inputExpression), error = function(x) {
		if(!any(trim(x$message) == trim(message))) {
			stop(paste("An error was thrown with the wrong message:",
				x$message), call. = FALSE)
		} else {
			assign("checkErrorState", TRUE, pkg_globals)
		}
	})
	if (!get("checkErrorState", pkg_globals)) {
		stop(paste("No error was observed for the expression",
			deparse(inputExpression, width.cutoff = 500L)), call. = FALSE)
	} else if (getOption("mxPrintUnitTests")) {
		cat(paste("The expected error was observed for the expression",
			deparse(inputExpression, width.cutoff = 500L), '\n'))
	}
}
