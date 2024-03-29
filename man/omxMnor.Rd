%
%   Copyright 2007-2018 by the individuals mentioned in the source code history
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

\name{omxMnor}
\alias{omxMnor}

\title{Multivariate Normal Integration}

\description{
   Given a covariance matrix, a means vector, and vectors of lower and upper bounds, returns the multivariate normal integral across the space between bounds.
}

\usage{
omxMnor(covariance, means, lbound, ubound)
}

\arguments{
   \item{covariance}{the covariance matrix describing the multivariate normal distribution.}
   \item{means}{a row vector containing means of the variables of the underlying distribution.}
   \item{lbound}{a row vector containing the lower bounds of the integration in each variable.}
   \item{ubound}{a row vector containing the upper bounds of the integration in each variable.}
}

\details{
   The order of columns in the \sQuote{means}, \sQuote{lbound}, and \sQuote{ubound} vectors are assumed to be the same as that of the covariance matrix.  That is, means[i] is considered to be the mean of the variable whose variance is in covariance[i,i].  That variable will be integrated from lbound[i] to ubound[i] as part of the integration.

   The value of ubound[i] or lbound[i] may be set to Inf or -Inf if a boundary at positive or negative infinity is desired.
   
   For all i, ubound[i] must be strictly greater than lbound[i].

   The algorithm for multivariate normal integration we use is Alan Genz's FORTRAN implementation of the SADMVN routine described by Genz (1992).
}

\references{
Genz, A.  (1992).  Numerical Computation of Multivariate Normal Probabilities.  \emph{Journal of Computational Graphical Statistics, 1,}  141-149.
}

\examples{

data(myFADataRaw)

covariance <- cov(myFADataRaw[,1:3])
means <- colMeans(myFADataRaw[,1:3])
lbound <- c(-Inf, 0,   1)    # Integrate from -Infinity to 0 on first variable 
ubound <- c(0,    Inf, 2.5)  # From 0 to +Infinity on second, and from 1 to 2.5 on third
omxMnor(covariance, means, lbound, ubound)
# 0.0005995

# An alternative specification of the bounds follows
# Integrate from -Infinity to 0 on first variable 
v1bound = c(-Inf, 0)
# From 0 to +Infinity on second
v2bound = c(0, Inf)
# and from 1 to 2.5 on third
v3bound = c(1, 2.5)
bounds <- cbind(v1bound, v2bound, v3bound)
lbound <- bounds[1,]  
ubound <- bounds[2,]  
omxMnor(covariance, means, lbound, ubound)

}
