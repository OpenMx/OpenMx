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

\name{rvectorize}
\alias{rvectorize}

\title{Vectorize By Row}

\description{
   This function returns the vectorization of an input matrix in a row by row traversal of the matrix.  The output is returned as a column vector.
}

\usage{
rvectorize(x)
}

\arguments{
   \item{x}{an input matrix.}
}

\seealso{
\code{\link{cvectorize}}, \code{\link{vech}}, \code{\link{vechs}} 
}

\examples{

rvectorize(matrix(1:9, 3, 3))
rvectorize(matrix(1:12, 3, 4))


}
