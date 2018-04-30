\name{mxComputeGenSA}
\alias{mxComputeGenSA}
\alias{MxComputeGenSA-class}
\title{Optimization using generalized simulated annealing}
\usage{
mxComputeGenSA(freeSet=NA_character_, ..., fitfunction='fitfunction',
 plan=mxComputeOnce('fitfunction','fit'),
 verbose=0L, qv=2.62, qaInit=-3., lambda=0.85,
 tempStart=5230., tempEnd=0.1, stepsPerTemp=1L)
}
\arguments{
\item{freeSet}{names of matrices containing free variables}

\item{...}{Not used.  Forces remaining arguments to be specified by name.}

\item{fitfunction}{name of the fitfunction (defaults to 'fitfunction')}

\item{plan}{compute plan to optimize the model}

\item{verbose}{level of debugging output}

\item{qv}{Controls the shape of the visiting distribution}

\item{qaInit}{Controls the shape of the initial acceptance distribution}

\item{lambda}{When nonzero, reduces the probability of acceptance in time}

\item{tempStart}{Starting temperature}

\item{tempEnd}{Ending temperature}

\item{stepsPerTemp}{Number of Markov steps per temperature step}
}
\description{
Performs generalized simulated annealing (GSA) to minimize the fit
function. The number of function evaluations are determined by the
\code{tempStart} and \code{tempEnd} parameters. There is no provision to
stop early because there is no way to determine whether the algorithm
has converged.
}
\details{
If the original starting values are outside of the feasible set,
a few attempts are made to find viable starting values.

Classical simulated annealing (CSA) can be obtained with
\code{qv=qa=1}. Fast simulated annealing (FSA) can be obtained with
\code{qv=2} and \code{qa=1}. FSA is faster than CSA, but GSA is
faster than FSA.
GenSA default parameters are set to those identified in 
Xiang, Sun, Fan & Gong (1997).  
}
\references{
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
	mxComputeGenSA())

m1 <- mxRun(m1)
summary(m1)
}
