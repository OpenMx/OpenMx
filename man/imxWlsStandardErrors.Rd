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

\name{imxWlsStandardErrors}
\alias{imxWlsStandardErrors}

\title{Calculate Standard Errors for a WLS Model}

\description{
   This is an internal function used to calculate standard errors for weighted least squares models.
}

\usage{
imxWlsStandardErrors(model)
}


\arguments{
   \item{model}{An MxModel object with acov (WLS) data}
}

\details{
The standard errors for models fit with maximum likelihood are related to the second derivative (Hessian) of the likelihood function with respect to the free parameters.  For models fit with weighted least squares a different expression is required.  If \eqn{J} is the first derivative (Jacobian) of the mapping from the free parameters to the unique elements of the expected covariance, means, and threholds, \eqn{V} is the weight matrix used, \eqn{W} is the inverse of the full weight matrix, and \eqn{U= V J (J' V J)^{-1}}, then the asymptotic covariance matrix of the free parameters is
\deqn{Acov(\theta) = U' W U}
with \eqn{U'} indicating the transpode of \eqn{U}.
}

\value{
A named list with components 
\describe{
	\item{SE}{The standard errors of the free parameters}
	\item{Cov}{The full covariance matrix of the free parameters.  The square root of the diagonal elements of Cov equals SE.}
}
}

\references{
M. W. Browne. (1984).  Asymptotically Distribution-Free Methods for the Analysis of Covariance Structures. \emph{British Journal of Mathematical and Statistical Psychology, 37}, 62-83.

F. Yang-Wallentin, K. G. J\ifelse{latex}{\out{\"o}}{\ifelse{html}{\out{&ouml;}}{o}}reskog, & H. Luo.  (2010).   Confirmatory Factor Analysis of Ordinal Variables with Misspecified Models.  \emph{Structural Equation Modeling, 17}, 392-423.
}


