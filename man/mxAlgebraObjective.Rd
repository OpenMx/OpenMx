%
%   Copyright 2007-2021 by the individuals mentioned in the source code history
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

\name{mxAlgebraObjective}
\alias{mxAlgebraObjective}

\title{DEPRECATED:  Create MxAlgebraObjective Object}

\description{
WARNING: Objective functions have been deprecated as of OpenMx 2.0.  

Please use MxFitFunctionAlgebra() instead.  As a temporary workaround, MxAlgebraObjective returns a list containing a NULL MxExpectation object and an MxFitFunctionAlgebra object.

All occurrences of

mxAlgebraObjective(algebra, numObs = NA, numStats = NA)

Should be changed to

mxFitFunctionAlgebra(algebra, numObs = NA, numStats = NA)
}

\arguments{
   \item{algebra}{A character string indicating the name of an \link{MxAlgebra} or \link{MxMatrix} object to use for optimization.}
   \item{numObs}{(optional) An adjustment to the total number of observations in the model.}
   \item{numStats}{(optional) An adjustment to the total number of observed statistics in the model.}
}

\details{
NOTE: THIS DESCRIPTION IS DEPRECATED.  Please change to using \link{mxFitFunctionAlgebra} as shown in the example below.

Fit functions are functions for which free parameter values are chosen such that the value of the objective function is minimized. While the other fit functions in OpenMx require an expectation function for the model, the \code{mxAlgebraObjective} function uses the referenced \code{\link{MxAlgebra}} or \code{\link{MxMatrix}} object as the function to be minimized.

If a model's primary objective function is a \code{mxAlgebraObjective} objective function, then the referenced algebra in the objective function must return a 1 x 1 matrix (when using OpenMx's default optimizer).  There is no restriction on the dimensions of an objective function that is not the primary, or \sQuote{topmost}, objective function.

To evaluate an algebra objective function, place the following objects in a \code{\link{MxModel}} object: a \code{MxAlgebraObjective}, \code{\link{MxAlgebra}} and \code{\link{MxMatrix}} entities referenced by the \code{MxAlgebraObjective}, and optional \code{\link{MxBounds}} and \code{\link{MxConstraint}} entities. This model may then be evaluated using the \code{\link{mxRun}} function. The results of the optimization may be obtained using the \code{\link{mxEval}} function on the name of the \code{\link{MxAlgebra}}, after the model has been run.
}

\value{
Returns a list containing a NULL MxExpectation object and an MxFitFunctionAlgebra object.  MxFitFunctionAlgebra objects should be included with models with referenced \code{\link{MxAlgebra}} and \code{\link{MxMatrix}} objects.
}

\references{
The OpenMx User's guide can be found at \url{https://openmx.ssri.psu.edu/documentation/}.
}

\seealso{
\link{mxAlgebra} to create an algebra suitable as a reference function to be minimized. More information about the OpenMx package may be found \link[=OpenMx]{here}. 
}

\examples{

# Create and fit a very simple model that adds two numbers using mxFitFunctionAlgebra

library(OpenMx)

# Create a matrix 'A' with no free parameters
A <- mxMatrix('Full', nrow = 1, ncol = 1, values = 1, name = 'A')

# Create an algebra 'B', which defines the expression A + A
B <- mxAlgebra(A + A, name = 'B')

# Define the objective function for algebra 'B'
objective <- mxFitFunctionAlgebra('B')

# Place the algebra, its associated matrix and 
# its objective function in a model
tmpModel <- mxModel(model="Addition", A, B, objective)

# Evalulate the algebra
tmpModelOut <- mxRun(tmpModel)

# View the results
tmpModelOut$output$minimum
}
