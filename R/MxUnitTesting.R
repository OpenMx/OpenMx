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
		stop(paste(omxQuotes(paste(a, collapse = ' ')), 
			"and", omxQuotes(paste(b, collapse = ' ')), 
			"are not both vectors")) 	
	}
	if((is.matrix(a) && (nrow(a) > 1 || ncol(a) > 1) && !is.matrix(b)) || 
		(is.matrix(b) && (nrow(b) > 1 || ncol(b) > 1) && !is.matrix(a))) {
		stop(paste(omxQuotes(paste(a, collapse = ' ')), 
			"and", omxQuotes(paste(b, collapse = ' ')), 
			"are not both matrices")) 	
	}	
	if (is.vector(a) && (length(a) != length(b))) {
		stop(paste(omxQuotes(paste(a, collapse = ' ')), 
			"and", omxQuotes(paste(b, collapse = ' ')), 
			"do not have equal length :",
			length(a), 'and', length(b)))
	}
	if (is.matrix(a) && (nrow(a) > 1 || ncol(a) > 1) && any(dim(a) != dim(b))) {
		stop(paste(omxQuotes(paste(a, collapse = ' ')), 
			"and", omxQuotes(paste(b, collapse = ' ')), 
			"are not of equal dimension :", 
			paste(dim(a), collapse = ' x '), 'and', 
			paste(dim(b), collapse = ' x ')))
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


omxCheckCloseEnough <- function(a, b, epsilon = 10^(-15)) {
	checkEqualDimensions(a, b)
	if(any(mapply(function(x,y) {
			is.na(a) || is.na(b) },
			as.vector(a), as.vector(b)))) {
		stop("omxCheckCloseEnough does not support NA values")
	}
	check <- any(mapply(function(x,y) {
			abs(x - y) > epsilon }, 
			as.vector(a), as.vector(b)))
	if (check) {
		stop(paste("In", deparse(width.cutoff = 400L, sys.call()), ":",
			omxQuotes(paste(a, collapse = ' ')), 
			"and", omxQuotes(paste(b, collapse = ' ')), 
			"are not equal to within", epsilon), 
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

omxCheckWarning <- function(expression, message) {
	inputExpression <- match.call()$expression
	checkWarningState <- FALSE
	tryCatch(eval(inputExpression), warning = function(x) {
		if(trim(x$message) != trim(message)) {
			stop(paste("An warning was thrown with the wrong message:",
				x$message), call. = FALSE)
		} else { checkWarningState <<- TRUE }
	})
	if (!checkWarningState) {
		stop(paste("No warning was observed for the expression",
			deparse(inputExpression, width.cutoff = 500L)), call. = FALSE)
	} else if (getOption("mxPrintUnitTests")) {
		cat(paste("The expected warning was observed for the expression",
			deparse(inputExpression, width.cutoff = 500L), '\n'))
	}
}

omxCheckError <- function(expression, message) {
	inputExpression <- match.call()$expression
	checkErrorState <- FALSE
	tryCatch(eval(inputExpression), error = function(x) {
		if(trim(x$message) != trim(message)) {
			stop(paste("An error was thrown with the wrong message:",
				x$message), call. = FALSE)
		} else { checkErrorState <<- TRUE }
	})
	if (!checkErrorState) {
		stop(paste("No error was observed for the expression",
			deparse(inputExpression, width.cutoff = 500L)), call. = FALSE)
	} else if (getOption("mxPrintUnitTests")) {
		cat(paste("The expected error was observed for the expression",
			deparse(inputExpression, width.cutoff = 500L), '\n'))
	}
}
