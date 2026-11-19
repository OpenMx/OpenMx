%
%   Copyright 2007-2026 by the individuals mentioned in the source code history
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

\name{mxComputeConfidenceInterval}
\alias{mxComputeConfidenceInterval}
\alias{MxComputeConfidenceInterval-class}
\title{Find likelihood-based confidence intervals}
\usage{
mxComputeConfidenceInterval(
  plan,
  ...,
  freeSet = NA_character_,
  verbose = 0L,
  engine = NULL,
  fitfunction = "fitfunction",
  tolerance = NA_real_,
  constraintType = "none"
)
}
\arguments{
\item{plan}{compute plan to optimize the model}

\item{...}{Not used.  Forces remaining arguments to be specified by name.}

\item{freeSet}{names of matrices containing free variables}

\item{verbose}{integer. Level of run-time diagnostic output. Set to zero to disable}

\item{engine}{\lifecycle{deprecated}}

\item{fitfunction}{the name of the deviance function}

\item{tolerance}{\lifecycle{deprecated}}

\item{constraintType}{one of c('ineq', 'none')}
}
\description{
OpenMx has two internal representations of the confidence-limit optimization problem: an inequality-\emph{constrained} representation and an \emph{unconstrained} (quadratic-penalty) representation.  Let \eqn{\theta} denote the vector of free parameters, \eqn{\hat{\theta}} denote the MLE, \eqn{L(\theta)} denote -2 times the loglikelihood evaluated at \eqn{\theta}, \eqn{\chi^2 _{crit}} denote the quantile of the 1\emph{df} chi-square corresponding to the desired coverage probability (e.g., 3.841 for 95\% coverage), and \eqn{\phi} denote a scalar-valued function of \eqn{\theta} for which a confidence limit is to be estimated.  Let \eqn{c=1} when searching for a lower confidence limit, and let \eqn{c=-1} when searching for an upper confidence limit.

The constrained optimization problem is\cr
minimize \eqn{c\phi}\cr
subject to \eqn{-L(\theta)+L(\hat{\theta})+{\chi^2 _{crit}}<0}\cr
.  The unconstrained optimization problem, which uses a quadratic penalty function to approximate an equality constraint, is\cr
minimize \eqn{\left(L(\theta)+{\chi^2 _{crit}}-L(\hat{\theta})\right)^2 + c\phi}\cr
.  It may be tempting to place a large coefficient in front of the quadratic penalty term to better approximate the equality constraint.  Unfortunately, the unconstrained optimization problem is usually ill-conditioned as it is, and introducing a penalty coefficient greater than unity exacerbates this ill-conditioning.\cr
The constrained representation is the theoretical superior of the two.  It lacks the biases of the unconstrained representation discussed in Neale & Miller (1997).  By default, SLSQP uses the constrained representation, whereas CSOLNP and NPSOL use the unconstrained representation.\cr
\subsection{Miscellaneous}{
\enumerate{
	\item{When estimating a confidence limit for a free parameter (and not an \link[=mxAlgebra]{algebra}), SLSQP automatically uses the (obvious) analytic gradient of the objective function, \eqn{c\phi}.}
	\item{SLSQP will re-attempt confidence-limit optimization up to three times, with start values randomly perturbed between attempts (\emph{a la} \code{\link{mxTryHard}()}).}
	\item{The Wu-Neale (2012) correction for a boundary on a free parameter only works correctly with the constrained representation.}
	}
}
}
\references{
Neale, M. C. & Miller M. B. (1997). The use of likelihood based
confidence intervals in genetic models.  \emph{Behavior Genetics,
27}(2), 113-120.

Pek, J. & Wu, H. (2015). Profile likelihood-based confidence intervals and regions for structural equation models.
\emph{Psychometrika, 80}(4), 1123-1145.

Wu, H. & Neale, M. C. (2012). Adjusted confidence intervals for a
bounded parameter. \emph{Behavior genetics, 42}(6), 886-898.
}
