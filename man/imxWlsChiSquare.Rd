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

\name{imxWlsChiSquare}
\alias{imxWlsChiSquare}

\title{Calculate Chi Square for a WLS Model}

\description{
   This is an internal function used to calculate the Chi Square distributed fit statistic for weighted least squares models.
}

\usage{
imxWlsChiSquare(model, J=NA)
}


\arguments{
   \item{model}{An MxModel object with acov (WLS) data}
   \item{J}{Optional pre-computed Jacobian matrix}
}

\details{
The Chi Square fit statistic for models fit with maximum likelihood depends on the difference in model fit in minus two log likelihood units between the saturated model and the more restricted model under investigation.  For models fit with weighted least squares a different expression is required.  If \eqn{J} is the first derivative (Jacobian) of the mapping from the free parameters to the unique elements of the expected covariance, means, and threholds, \eqn{J_c} is the orthogonal complement of \eqn{J}, \eqn{W} is the inverse of the full weight matrix, and \eqn{e} is the difference between the sample-estimated and model-implied covariance, means, and thresholds, then the Chi Square fit statistic is
\deqn{\chi^2 = e' J_c (J'_c W J_c)^-1 J'_c e}
with \eqn{e'} indicating the transpose of \eqn{e}.  This Equation 2.20a from Browne (1984) where he showed that this statistic is chi-square distributed with the conventional degrees of freedom.

Mean and variance adjusted Chi Square statistics are also computed following Asparouhov and Muthen (2006).
}

\value{A named list with components
\describe{
	\item{Chi}{numeric value of the Chi Square fit statistic.}
	\item{ChiDoF}{degrees of freedom for the Chi Square fit statistic.}
	\item{ChiM}{numeric value of the mean adjusted Chi Square fit statistic}
	\item{ChiMV}{numeric value of the mean and variance adjusted Chi Square fit statistic}
	\item{mAdjust}{numeric value of the mean adjustment}
	\item{mvAdjust}{numeric value of the mean and variance adjustment}
	\item{dstar}{adjusted degrees of freedom for the mean and variance adjusted Chi Square fit statistic}
}
}

\references{
M. W. Browne. (1984).  Asymptotically Distribution-Free Methods for the Analysis of Covariance Structures. \emph{British Journal of Mathematical and Statistical Psychology}, \strong{37}, 62-83.

T. Asparouhov and B. O. Muthen.  (2006).  Robust Chi Square Difference Testing with Mean and Variance Adjusted Test Statistics.  \emph{Mplus Web Notes: No. 10}.
}

