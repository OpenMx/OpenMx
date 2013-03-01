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

\name{vechs}
\alias{vechs}

\title{Strict Half-vectorization}

\description{
   This function returns the strict half-vectorization of an input matrix as a column vector.
}

\usage{
vechs(x)
}

\arguments{
   \item{x}{an input matrix.}
}

\details{
   The half-vectorization of an input matrix consists of the elements in the lower triangle of the matrix, excluding the elements along the diagonal of the matrix, as a column vector.  The column vector is created by traversing the matrix in column-major order.
}

\seealso{
\code{\link{vech}}, \code{\link{rvectorize}}, \code{\link{cvectorize}}
}

\examples{

vechs(matrix(1:9, 3, 3))
vechs(matrix(1:12, 3, 4))


}
