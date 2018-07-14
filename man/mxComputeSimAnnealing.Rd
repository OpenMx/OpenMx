\name{mxComputeSimAnnealing}
\alias{mxComputeSimAnnealing}
\alias{MxComputeSimAnnealing-class}
\title{Optimization using generalized simulated annealing}
\usage{
mxComputeSimAnnealing(freeSet=NA_character_, ..., fitfunction='fitfunction',
 plan=mxComputeOnce('fitfunction','fit'),
 verbose=0L, method=c("tsallis1996", "ingber2012"), control=list(),
 defaultGradientStepSize=imxAutoOptionValue("Gradient step size"),
 defaultFunctionPrecision=imxAutoOptionValue("Function precision"))
}
\arguments{
\item{freeSet}{names of matrices containing free variables}

\item{...}{Not used.  Forces remaining arguments to be specified by name.}

\item{fitfunction}{name of the fitfunction (defaults to 'fitfunction')}

\item{plan}{compute plan to optimize the model}

\item{verbose}{level of debugging output}

\item{method}{which algorithm to use}

\item{control}{control parameters specific to the chosen method}

\item{defaultGradientStepSize}{the default gradient step size}

\item{defaultFunctionPrecision}{the default function precision}
}
\description{
Performs simulated annealing to minimize the fit function.
If the original starting values are outside of the feasible set,
a few attempts are made to find viable starting values.
}
\details{
For method \sQuote{tsallis1996},
the number of function evaluations are determined by the
\code{tempStart} and \code{tempEnd} parameters. There is no provision to
stop early because there is no way to determine whether the algorithm
has converged. The Markov step is implemented by cycling through each
parameters in turn and considering a univariate jump (like a Gibbs sampler).

Control parameters include \code{qv} to control the shape of the
visiting distribution, \code{qaInit} to control the shape of the initial
acceptance distribution, \code{lambda} to reduce the probability of
acceptance in time, \code{tempStart} to specify starting temperature,
\code{tempEnd} to specify ending temperature, and \code{stepsPerTemp} to
set the number of Markov steps per temperature step.

Non-linear constraints are accommodated by a penalty function.
Inequality constraints work reasonably well, but
equality constraints do not work very well.
Constrained optimization will likely require increasing \code{stepsPerTemp}.

Classical simulated annealing (CSA) can be obtained with
\code{qv=qa=1} and \code{lambda=0}.
Fast simulated annealing (FSA) can be obtained with
\code{qv=2}, \code{qa=1}, and \code{lambda=0}.
FSA is faster than CSA, but GSA is faster than FSA.
GenSA default parameters are set to those identified in 
Xiang, Sun, Fan & Gong (1997).  

Method \sQuote{ingber2012} has spawned a cultural tradition over more
than 30 years that is documented in Aguiar e Oliveira et al (2012).
Options are specified using the traditional option names in the
\code{control} list. However, there are a few option changes to
make ASA fit better with OpenMx.
Instead of option \code{Curvature_0}, use \link{mxComputeNumericDeriv}.
ASA_PRINT output is directed to \code{/dev/null} by default.
To direct ASA_PRINT output to console use \code{control=list('Asa_Out_File'= '/dev/fd/1')}.
ASA's option to control the finite differences gradient step size,
\code{Delta_X}, defaults to \link{mxOption}'s \sQuote{Gradient step
  size} instead of ASA's traditional 0.001.
Similarly, \code{Cost_Precision} defaults to \link{mxOption}'s
\sQuote{Function Precision} instead of ASA's traditional 1e-18.
}
\references{
Aguiar e Oliveira, H., Ingber, L., Petraglia, A., Petraglia, M. R., & Machado, M. A. S. (2012). \emph{Stochastic global optimization and its applications with fuzzy adaptive simulated annealing.} Springer Publishing Company, Incorporated.

Tsallis, C., & Stariolo, D. A. (1996). Generalized simulated
annealing. \emph{Physica A: Statistical Mechanics and its Applications,
233}(1-2), 395-406.

Xiang, Y., Sun, D. Y., Fan, W., & Gong, X. G. (1997). Generalized
simulated annealing algorithm and its application to the Thomson
model. \emph{Physics Letters A, 233}(3), 216-220.
}
\seealso{
\code{\link{mxComputeTryHard}}
}
\examples{
library(OpenMx)
m1 <- mxModel(
	"poly22", # Eqn 22 from Tsallis & Stariolo (1996)
	mxMatrix(type='Full', values=runif(4, min=-1e6, max=1e6),
		ncol=1, nrow=4, free=TRUE, name='x'),
	mxAlgebra(sum((x*x-8)^2) + 5*sum(x) + 57.3276, name="fit"),
	mxFitFunctionAlgebra('fit'),
	mxComputeSimAnnealing())

m1 <- mxRun(m1)
summary(m1)
}
