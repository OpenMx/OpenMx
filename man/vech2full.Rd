%
%   Copyright 2007-2015 The OpenMx Project
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

\name{vech2full}
\alias{vech2full}

\title{Inverse Half-vectorization}

\description{
   This function returns the symmetric matrix constructed from a half-vectorization.
}

\usage{
vech2full(x)
}

\arguments{
   \item{x}{an input single column or single row matrix.}
}

\details{
   The half-vectorization of an input matrix consists of the elements in the lower triangle of the matrix, including the elements along the diagonal of the matrix, as a column vector.  The column vector is created by traversing the matrix in column-major order.  The inverse half-vectorization takes a vector and reconstructs a symmetric matrix such that \code{vech2full(vech(x))} is identical to \code{x} if \code{x} is symmetric.

   Note that very few vectors have the correct number of elements to construct a symmetric matrix.  For example, vectors with 1, 3, 6, 10, and 15 elements can be used to make a symmetric matrix, but none of the other numbers between 1 and 15 can.  An error is thrown if the number of elements in \code{x} cannot be used to make a symmetric matrix.
}

\seealso{
\code{\link{vechs2full}}, \code{\link{vech}}, \code{\link{vechs}}, \code{\link{rvectorize}}, \code{\link{cvectorize}}
}

\examples{

vech2full(1:10)

matrix(1:16, 4, 4)
vech(matrix(1:16, 4, 4))
vech2full(vech(matrix(1:16, 4, 4)))



}
