%
%   Copyright 2007-2014 The OpenMx Project
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

\name{mxFactor}
\alias{mxFactor}

\title{Fail-safe Factors}

\description{
   This is a wrapper for the R function \code{\link{factor}}. 

OpenMx requires ordinal data to be ordered. R's \code{factor} function doesn't enforce this, hence this wrapper exists to throw an error should you accidentally try and run with ordered = FALSE.

Also, the \sQuote{levels} parameter is optional in R's \code{factor} function. However, relying on the data to specify the data is foolhardy for the following reasons: The \code{factor} function will skip levels missing from the data: Specifying these in levels leaves the list of levels complete. Data will often not explore the min and max level that the user knows are possible. For these reasons this function forces you to write out all possible levels explicitly.

}

\usage{
mxFactor(x = character(), levels, labels = levels, 
    exclude = NA, ordered = TRUE)
}

\arguments{
    \item{x}{either a vector of data or a data.frame object.}
    \item{levels}{a mandatory vector of the values that 'x' might have taken.}
    \item{labels}{_either_ an optional vector of labels for the levels, _or_ a character string of length 1.}
    \item{exclude}{a vector of values to be excluded from the set of levels.}
    \item{ordered}{logical flag to determine if the levels should be regarded as ordered (in the order given).  Required to be TRUE.}
}

\details{
If \sQuote{x} is a data.frame, then all of the columns of \sQuote{x} are converted into ordered factors. If \sQuote{x} is a data.frame, then \sQuote{levels} and \sQuote{labels} may be either a list or a vector. When \sQuote{levels} is a list, then different levels are assigned to different columns of the constructed data.frame object. When \sQuote{levels} is a vector, then the same levels are assigned to all the columns of the data.frame object. The function will throw an error if \sQuote{ordered} is not TRUE or if \sQuote{levels} is missing.  See \code{\link{factor}} for more information on creating ordered factors. 

}

\references{
The OpenMx User's guide can be found at http://openmx.psyc.virginia.edu/documentation.
}

\examples{
myVar <- c("s", "t", "a", "t", "i", "s", "t", "i", "c", "s")
ff    <- mxFactor(myVar, levels=letters)  # letters is a built in list of all lowercase letters of the alphabet
ff
# [1] s t a t i s t i c s
# Levels: a < b < c < d < e < f < g < h < i < j < k < l < m < n < o < p < q < r < s < t < u < v < w < x < y < z

as.integer(ff)  # the internal codes

factor(ff)      # NOTE: drops the levels that do not occur.
                # mxFactor prevents you doing this unintentionally.

# This example works on a dataframe
foo <- data.frame(x=c(1:3),y=c(4:6),z=c(7:9))
mxFactor(foo, c(1:9))                       # Applys one set of levels to all three columns
mxFactor(foo, list(c(1:3), c(4:6), c(7:9))) # Apply unique sets of levels to each variable
    
}
