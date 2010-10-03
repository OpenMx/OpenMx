%
%   Copyright 2007-2010 The OpenMx Project
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

\name{vec2diag}
\alias{vec2diag}

\title{Create Diagonal Matrix From Vector}

\description{
   Given an input row or column vector, \code{vec2diag} returns a diagonal matrix with the input argument along the diagonal.
}

\usage{
vec2diag(x)
}

\arguments{
   \item{x}{a row or column vector.}
}

\details{
   Similar to the function \code{\link{diag}}, except that the input argument is always
   treated as a vector of elements to place along the diagonal.
}

\seealso{
\code{\link{diag2vec}} 
}

\examples{

vec2diag(matrix(1:4, 1, 4))
vec2diag(matrix(1:4, 4, 1))


}
