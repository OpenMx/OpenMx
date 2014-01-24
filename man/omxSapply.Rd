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

\name{omxSapply}
\alias{omxSapply}

\title{On-Demand Parallel Sapply}

\description{
   If the snowfall library is loaded, then this function calls \code{\link[snowfall]{sfSapply}}.  Otherwise it invokes \code{\link{sapply}}.
}

\usage{
omxSapply(x, fun, ..., simplify = TRUE, USE.NAMES = TRUE)
}

\arguments{
   \item{x}{a vector (atomic or list) or an expressions vector. Other objects (including classed objects) will be coerced by \code{\link{as.list}}.}
   \item{fun}{the function to be applied to each element of \code{x}.}
   \item{...}{optional arguments to \code{fun}.}
   \item{simplify}{logical; should the result be simplified to a vector or matrix if possible?}
   \item{USE.NAMES}{logical; if \code{TRUE} and if \code{x} is a character, use \code{x} as \code{\link{names}} for the result unless it had names already.}
}

\seealso{
\code{\link{omxApply}}, \code{\link{omxLapply}} 
}

\examples{

x <- list(a = 1:10, beta = exp(-3:3), logic = c(TRUE,FALSE,FALSE,TRUE))
# compute the list mean for each list element
omxSapply(x, quantile)

}
