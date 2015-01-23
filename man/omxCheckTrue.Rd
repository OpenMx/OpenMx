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

\name{omxCheckTrue}
\alias{omxCheckTrue}

\title{Boolean Equality Testing Function}

\description{
   This function tests whether an object is equal to TRUE.
}

\usage{
omxCheckTrue(a)
}

\arguments{
   \item{a}{the value to test.}
}

\details{
Checks element-wise whether an object is equal to TRUE. If any of the elements are false, then an error will be thrown. If \sQuote{a} is TRUE, by default the function will print a statement informing the user the test has passed.  To turn off these print statements use \code{options("mxPrintUnitTests" = FALSE)}.
}

\references{
The OpenMx User's guide can be found at http://openmx.psyc.virginia.edu/documentation.
}

\seealso{
\code{\link{omxCheckCloseEnough}}, \code{\link{omxCheckWithinPercentError}}, \code{\link{omxCheckIdentical}}, \code{\link{omxCheckSetEquals}}, \code{\link{omxCheckEquals}}
}

\examples{

omxCheckTrue(1 + 1 == 2)

omxCheckTrue(matrix(TRUE, 3, 3))

# Throws an error
try(omxCheckTrue(FALSE))

}
