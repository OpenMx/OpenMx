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

\name{omxCheckEquals}
\alias{omxCheckEquals}

\title{Equality Testing Function}

\description{
   This function tests whether two objects are equal using the \sQuote{==} operator.
}

\usage{
omxCheckEquals(...)
}

\arguments{
   \item{...}{arguments forwarded to \link[testthat:equality-expectations]{expect_equivalent}}
}

\details{
Performs the \sQuote{==} comparison on the two arguments. If the two arguments are not equal, then an error will be thrown.
}

\references{
The OpenMx User's guide can be found at http://openmx.ssri.psu.edu/documentation.
}

\seealso{
\code{\link{omxCheckCloseEnough}}, \code{\link{omxCheckWithinPercentError}}, \code{\link{omxCheckSetEquals}}, \code{\link{omxCheckTrue}}, \code{\link{omxCheckIdentical}}
}

\examples{

omxCheckEquals(c(1, 2, 3), c(1, 2, 3))

omxCheckEquals(FALSE, FALSE)

# Throws an error
try(omxCheckEquals(c(1, 2, 3), c(2, 1, 3)))

}
