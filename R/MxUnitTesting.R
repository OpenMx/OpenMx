#
#   Copyright 2007-2018 by the individuals mentioned in the source code history
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


omxCheckIdentical <- function(...) {
  if (requireNamespace('testthat', quietly = TRUE)) {
    testthat::expect_identical(...)
  } else {
    stop(paste0("Please install.packages(testthat) and try again"))
  }
}

omxCheckEquals <- function(...) {
  if (requireNamespace('testthat', quietly = TRUE)) {
    testthat::expect_equivalent(...)
  } else {
    stop(paste0("Please install.packages(testthat) and try again"))
  }
}

omxCheckSetEquals <- function(...) {
  if (requireNamespace('testthat', quietly = TRUE)) {
    testthat::expect_setequal(...)
  } else {
    stop(paste0("Please install.packages(testthat) and try again"))
  }
}

omxCheckTrue <- function(object) {
  if (requireNamespace('testthat', quietly = TRUE)) {
    if (is.numeric(object)) object <- as.logical(object)
    if (length(object) > 1) object <- all(object)
    testthat::expect_true(object)
  } else {
    stop(paste0("Please install.packages(testthat) and try again"))
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
##' be thrown.
##'
##' @param a a numeric vector or matrix
##' @param b a numeric vector or matrix
##' @param epsilon a non-negative tolerance threshold
##' @seealso
##' \code{\link{omxCheckWithinPercentError}},
##' \code{\link{omxCheckIdentical}}, \code{\link{omxCheckSetEquals}},
##' \code{\link{omxCheckTrue}}, \code{\link{omxCheckEquals}}
##' @references
##' The OpenMx User's guide can be found at http://openmx.ssri.psu.edu/documentation.
##' @examples
##' omxCheckCloseEnough(c(1, 2, 3), c(1.1, 1.9 ,3.0), epsilon = 0.5)
##' omxCheckCloseEnough(matrix(3, 3, 3), matrix(4, 3, 3), epsilon = 2)
##' # Throws an error
##' try(omxCheckCloseEnough(c(1, 2, 3), c(1.1, 1.9 ,3.0), epsilon = 0.01))
omxCheckCloseEnough <- function(a, b, epsilon = 10^(-15)) {
  if (requireNamespace('testthat', quietly = TRUE)) {
    testthat::expect_equivalent(object=a, expected=b, scale=1, tolerance=epsilon)
  } else {
    stop(paste0("Please install.packages(testthat) and try again"))
  }
}

omxCheckWithinPercentError <- function(a, b, percent = 0.1) {
  if (requireNamespace('testthat', quietly = TRUE)) {
    mapply(function(a1,b1) {
      testthat::expect_equivalent(object=a1, expected=b1, scale=abs(a1), tolerance=percent)
    }, as.vector(a), as.vector(b))
  } else {
    stop(paste0("Please install.packages(testthat) and try again"))
  }
}

##' Correct Warning Message Function
##'
##' This function tests whether the correct warning message is thrown.
##' Arguments \sQuote{expression} and \sQuote{message} give the expression
##' that generates the warning and the message that is supposed to be generated, respectively.
##' 
##' \emph{note}: to test for no warning, set \code{message = NA}.
##' 
##'
##' @param expression an R expression that produces a warning
##' @param message a character string with the desired warning message
##' @seealso
##' \code{\link{omxCheckError}}
##' \code{\link{omxCheckWithinPercentError}},
##' \code{\link{omxCheckIdentical}}, \code{\link{omxCheckSetEquals}},
##' \code{\link{omxCheckTrue}}, \code{\link{omxCheckEquals}}
##' @references
##' The OpenMx User's guide can be found at http://openmx.ssri.psu.edu/documentation.
##' @examples
##' foo <- omxCheckWarning(mxFIMLObjective('cov', 'mean'), "deprecated")
##' 
##' # Test for no warning
##' omxCheckWarning(2+2, message = NA)
##' 
omxCheckWarning <- function(expression, message) {
  if (requireNamespace('testthat', quietly = TRUE)) {
    testthat::expect_warning(expression, message, fixed=TRUE)
  } else {
    stop(paste0("Please install.packages(testthat) and try again"))
  }
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
##' The OpenMx User's guide can be found at http://openmx.ssri.psu.edu/documentation.
##' @examples
##' A <- mxMatrix('Full', 1, 1, labels = 'data.foo', free = TRUE, name = 'A')
##' model <- mxModel('model', A)
##' omxCheckError(mxRun(model), 
##'	paste("The definition variable 'data.foo'",
##'		"has been assigned to a",
##'		"free parameter in matrix 'A'"))
##' omxCheckCloseEnough(matrix(3, 3, 3), matrix(4, 3, 3), epsilon = 2)
omxCheckError <- function(expression, message) {
  if (requireNamespace('testthat', quietly = TRUE)) {
    testthat::expect_error(expression, message, fixed=TRUE)
  } else {
    stop(paste0("Please install.packages(testthat) and try again"))
  }
}
