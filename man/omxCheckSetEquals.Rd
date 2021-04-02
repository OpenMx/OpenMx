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

\name{omxCheckSetEquals}
\alias{omxCheckSetEquals}

\title{Set Equality Testing Function}

\description{
   This function tests whether two vectors contain the same elements.
}

\usage{
omxCheckSetEquals(...)
}

\arguments{
   \item{...}{arguments forwarded to \link[testthat:equality-expectations]{expect_setequal}}
}

\details{
Performs the \sQuote{setequal} function on the two arguments. If the two arguments do not contain the same elements, then an error will be thrown.
}

\references{
The OpenMx User's guide can be found at \url{https://openmx.ssri.psu.edu/documentation/}
}

\seealso{
\code{\link{omxCheckCloseEnough}}, \code{\link{omxCheckWithinPercentError}}, \code{\link{omxCheckIdentical}}, \code{\link{omxCheckTrue}}, \code{\link{omxCheckEquals}}
}

\examples{

omxCheckSetEquals(c(1, 1, 2, 2, 3), c(3, 2, 1))

omxCheckSetEquals(matrix(1, 1, 1), matrix(1, 3, 3))

# Throws an error
try(omxCheckSetEquals(c(1, 2, 3, 4), c(2, 1, 3)))

}
