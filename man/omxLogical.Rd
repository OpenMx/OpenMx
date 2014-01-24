%
%   Copyright 2007-2012 The OpenMx Project
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

\name{omxLogical}
\alias{omxLogical}
\alias{omxNot}
\alias{omxAnd}
\alias{omxOr}
\alias{omxGreaterThan}
\alias{omxLessThan}
\alias{omxApproxEquals}


\title{Logical mxAlgebra() operators}

\description{
   \code{omxNot} computes the unary negation of the values of a matrix.
   \code{omxAnd} computes the binary and of two matrices.
   \code{omxOr} computes the binary or of two matrices.
   \code{omxGreaterThan} computes a binary greater than of two matrices.
   \code{omxLessThan} computes the binary less than of two matrices.
   \code{omxApproxEquals} computes a binary equals within a specified epsilon of two matrices.
}

\usage{
omxNot(x)
omxAnd(x, y)
omxOr(x, y)
omxGreaterThan(x, y)
omxLessThan(x, y)
omxApproxEquals(x, y, epsilon)
}

\arguments{
   \item{x}{the first argument, the matrix which the logical operation will be applied to.}
   \item{y}{the second argument, applicable to binary functions.}
   \item{epsilon}{the third argument, specifies the error threshold for omxApproxEquals. Abs(x[i][j]-y[i][j]) must be less than epsilon[i][j].}
}

\examples{
A <- mxMatrix(values = runif(25), nrow = 5, ncol = 5, name = 'A')
B <- mxMatrix(values = runif(25), nrow = 5, ncol = 5, name = 'B')
EPSILON <- mxMatrix(values = 0.04*1:25, nrow = 5, ncol = 5, name = "EPSILON")

model <- mxModel(A, B, EPSILON, name = 'model')

mxEval(omxNot(A), model)
mxEval(omxGreaterThan(A,B), model)
mxEval(omxLessThan(B,A), model)
mxEval(omxOr(omxNot(A),B), model)
mxEval(omxAnd(omxNot(B), A), model)
mxEval(omxApproxEquals(A,B, EPSILON), model)
}
