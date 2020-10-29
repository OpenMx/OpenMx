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


omxCheckIdentical <- function(...) {
  if (requireNamespace('testthat', quietly = TRUE)) {
    testthat::expect_identical(...)
  } else {
    stop(paste0("Please install.packages(testthat) and try again")) #nocov
  }
}

omxCheckEquals <- function(...) {
  if (requireNamespace('testthat', quietly = TRUE)) {
    testthat::expect_equivalent(...)
  } else {
    stop(paste0("Please install.packages(testthat) and try again")) #nocov
  }
}

omxCheckSetEquals <- function(...) {
  if (requireNamespace('testthat', quietly = TRUE)) {
    testthat::expect_setequal(...)
  } else {
    stop(paste0("Please install.packages(testthat) and try again")) #nocov
  }
}

omxCheckTrue <- function(a) {
  if (requireNamespace('testthat', quietly = TRUE)) {
    if (is.numeric(a)) a <- as.logical(a) #nocov
    if (length(a) > 1) a <- all(a)
    testthat::expect_true(a)
  } else {
    stop(paste0("Please install.packages(testthat) and try again")) #nocov
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
    if (is(a, "logLik")) a <- as.numeric(a) #nocov
    if (is(b, "logLik")) b <- as.numeric(b) #nocov
    testthat::expect_equivalent(object=a, expected=b, scale=1, tolerance=epsilon)
  } else {
    stop(paste0("Please install.packages(testthat) and try again")) #nocov
  }
}

omxCheckWithinPercentError <- function(a, b, percent = 0.1) {
  if (requireNamespace('testthat', quietly = TRUE)) {
    mapply(function(a1,b1) {
      if (is(a1, "logLik")) a1 <- as.numeric(a1) #nocov
      if (is(b1, "logLik")) b1 <- as.numeric(b1) #nocov
      testthat::expect_equivalent(object=a1, expected=b1, scale=abs(a1), tolerance=percent)
    }, as.vector(a), as.vector(b))
  } else {
    stop(paste0("Please install.packages(testthat) and try again")) #nocov
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
    stop(paste0("Please install.packages(testthat) and try again")) #nocov
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
omxCheckError <- function(expression, message) {
  if (requireNamespace('testthat', quietly = TRUE)) {
    testthat::expect_error(expression, message, fixed=TRUE)
  } else {
    stop(paste0("Please install.packages(testthat) and try again")) #nocov
  }
}
