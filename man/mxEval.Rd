%
%   Copyright 2007-2013 The OpenMx Project
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

\name{mxEval}
\alias{mxEval}

\title{Evaluate Values in MxModel}

\description{
   This function can be used to evaluate an arbitrary R expression that includes named entities from a \link{MxModel} object, or labels from a \link{MxMatrix} object.
}

\usage{
mxEval(expression, model, compute = FALSE, show = FALSE, defvar.row = 1,
    cache = new.env(parent = emptyenv()), cacheBack = FALSE)
}

\arguments{
   \item{expression}{An arbitrary R expression.}
   \item{model}{The model in which to evaluate the expression.}
   \item{compute}{If TRUE then compute the value of algebra expressions.}
   \item{show}{If TRUE then print the translated expression.}
   \item{defvar.row}{The row number for definition variables when compute=TRUE.}
   \item{cache}{An R environment of matrix values used to speedup computation.}
   \item{cacheBack}{If TRUE then return the list pair (value, cache).}
}

\details{
The argument \sQuote{expression} is an arbitrary R expression.  Any named entities that are used within the R expression are translated into their current value from the model. Any labels from the matrices within the model are translated into their current value from the model. Finally the expression is evaluated and the result is returned.  To enable debugging, the \sQuote{show} argument has been provided.  The most common mistake when using this function is to include named entities in the model that are identical to R function names.  For example, if a model contains a named entity named \sQuote{c}, then the following mxEval call will return an error: \code{mxEval(c(A, B, C), model)}.

If \sQuote{compute} is FALSE, then MxAlgebra expressions return their current values as they have been computed by the optimization call (using \link{mxRun}). If the \sQuote{compute} argument is TRUE, then MxAlgebra expressions will be calculated in R.  Any references to an objective function that has not yet been calculated will return a 1 x 1 matrix with a value of NA.

The \sQuote{cache} is used to speedup calculation by storing previously computing values. The cache is a list of matrices, such that names(cache) must all be of the form \dQuote{modelname.entityname}. Setting \sQuote{cacheBack} 
to TRUE will return the pair list(value, cache) where value is the result of the mxEval() computation and cache is the updated cache.}

\references{
The OpenMx User's guide can be found at http://openmx.psyc.virginia.edu/documentation.
}

\seealso{
\link{mxAlgebra} to create algebraic expressions inside your model 
and \link{mxModel} for the model object mxEval looks inside when evaluating. 
}

\examples{

matrixA <- mxMatrix("Full", nrow = 1, ncol = 1, 
    values = 1, name = "A")
algebraB <- mxAlgebra(A + A, name = "B")

model <- mxModel(matrixA, algebraB, name = "test")
model <- mxRun(model)
results <- mxEval(B, model, compute=TRUE, cacheBack=TRUE)
B <- get("test.B", envir = results[[2]])
A <- get("test.A", envir = results[[2]])

\dontrun{

start <- mxEval(-pi * A, model) 
mxEval(plot(sin, start, B * pi), model)

# The statement above is equivalent to:

plot(sin, -pi, 2 * pi)

}

}
