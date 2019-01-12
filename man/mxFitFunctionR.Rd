%
%   Copyright 2007-2019 by the individuals mentioned in the source code history
%
%   Licensed under the Apache License, Version 2.0 (the "License");
%   you may not use this file except in compliance with the License.
%   You may obtain a copy of the License at
% 
%        http://www.apache.org/licenses/LICENSE-2.0
% 
%   Unless required by applicable law or agreed to in writing, software
%   distributed under the License is distributed on an "AS IS" BASIS,
%   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
%   See the License for the specific language governing permissions and
%   limitations under the License.

\name{mxFitFunctionR}
\alias{mxFitFunctionR}
\alias{MxFitFunctionR-class}
\alias{print,MxFitFunctionR-method}
\alias{show,MxFitFunctionR-method}

\title{Create MxFitFunctionR Object}

\description{
mxFitFunctionR returns an MxFitFunctionR object.
}

\usage{
mxFitFunctionR(fitfun, ..., units="-2lnL")
}

\arguments{
   \item{fitfun}{A function that accepts two arguments.}
   \item{...}{The initial state information to the objective function.}
   \item{units}{(optional) The units of the fit statistic.}
}

\details{
The mxFitFunctionR function evaluates a user-defined R function called the 'fitfun'.  mxFitFunctionR is useful in defining new mxFitFunctions, since any calculation that can be performed in R can be treated as an mxFitFunction.

The 'fitfun' argument must be a function that accepts two arguments. The first argument
is the mxModel that should be evaluated, and the second argument is some persistent 
state information that can be stored between one iteration of optimization to the next
iteration. It is valid for the function to simply ignore the second argument.

The function must return either a single numeric value, or a list of exactly two elements.
If the function returns a list, the first argument must be a single numeric value and the 
second element will be the new persistent state information to be passed into this function
at the next iteration. The single numeric value will be used by the optimizer to perform
optimization.

The initial default value for the persistent state information is NA.

Throwing an exception (via stop) from inside fitfun may result
in unpredictable behavior. You may want to wrap your code in
tryCatch while experimenting.

fitfun should not call R functions that use OpenMx's compiled backend, including (but not limited to) \code{\link{omxMnor}()}, because doing so can crash R.
}

\value{
Returns an MxFitFunctionR object. 
}

\seealso{
Other fit functions:
\code{\link{mxFitFunctionMultigroup}}, \code{\link{mxFitFunctionML}}, 
\code{\link{mxFitFunctionWLS}}, \code{\link{mxFitFunctionAlgebra}}, 
\code{\link{mxFitFunctionGREML}}, 
\code{\link{mxFitFunctionRow}}

More information about the OpenMx package may be found \link[=OpenMx]{here}. 
}

\references{
The OpenMx User's guide can be found at http://openmx.ssri.psu.edu/documentation.
}

\examples{

# Create and fit a model using mxFitFunctionR

library(OpenMx)

A <- mxMatrix(nrow = 2, ncol = 2, values = c(1:4), free = TRUE, name = 'A')
squared <- function(x) { x ^ 2 }

# Define the objective function in R

objFunction <- function(model, state) {
    values <- model$A$values 
    return(squared(values[1,1] - 4) + squared(values[1,2] - 3) +
        squared(values[2,1] - 2) + squared(values[2,2] - 1))
}

# Define the expectation function

fitFunction <- mxFitFunctionR(objFunction)

# Define the model

tmpModel <- mxModel(model="exampleModel", A, fitFunction)

# Fit the model and print a summary

tmpModelOut <- mxRun(tmpModel)
summary(tmpModelOut)

}
