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

\name{tr}
\alias{tr}

\title{trace}

\description{
   This function returns the trace of an n-by-n square matrix x,  defined to be the sum of the elements
   on the main diagonal (the diagonal from the upper left to the lower right).

}

\usage{
tr(x)
}

\arguments{
   \item{x}{an input matrix. Must be square}
}

\details{
The input matrix must be square.
}

\seealso{
\code{\link{vech}}, \code{\link{rvectorize}}, \code{\link{cvectorize}}
}

\examples{

tr(matrix(1:9, 3, 3))
tr(matrix(1:12, 3, 4))


}
