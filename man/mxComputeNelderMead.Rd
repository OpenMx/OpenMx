\name{mxComputeNelderMead}
\alias{mxComputeNelderMead}
\alias{MxComputeNelderMead}
\alias{MxComputeNelderMead-class}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Optimize parameters using a variation of the Nelder-Mead algorithm.
}
\description{
OpenMx includes a flexible, options-rich implementation of the Nelder-Mead algorithm.
}
\usage{
mxComputeNelderMead(
	freeSet=NA_character_, fitfunction="fitfunction", verbose=0L, 
	nudgeZeroStarts=mxOption(NULL,"Nudge zero starts"), 
	maxIter=NULL,	...,
	alpha=1, betao=0.5, betai=0.5, gamma=2, sigma=0.5, bignum=1e35, 
	iniSimplexType=c("regular","right","smartRight","random"),
	iniSimplexEdge=1, iniSimplexMat=NA, greedyMinimize=FALSE, 
	altContraction=FALSE, degenLimit=0, stagnCtrl=c(-1L,-1L),
	xTolProx=1e-4, fTolProx=1e-4,
	ineqConstraintMthd=c("soft","eqMthd"), 
	eqConstraintMthd=c("soft"))
}
\arguments{
  \item{freeSet}{Character-string names of \link[MxMatrix]{MxMatrices} containing free parameters.}
  \item{fitfunction}{Character-string name of the fitfunction; defaults to 'fitfunction'.}
  \item{verbose}{Integer level of reporting printed to terminal at \link[mxRun]{runtime}; defaults to 0.}
  \item{nudgeZeroStarts}{Should free parameters with start values of zero be "nudged" to 0.1 at \link[mxRun]{runtime}?  Defaults to the current global value of \link{mxOption} "Nudge zero starts".  May be a logical value, or one of character strings "Yes" or "No".}
  \item{maxIter}{Integer maximum number of iterations.  Value of \code{NULL} is accepted, in which case the value used at \link[mxRun]{runtime} will be 5 times the number of iterations specified by the effective value of \link{mxOption} "Major iterations".}
  \item{...}{Not used.  Forces remaining arguments to be specified by name.}
  \item{alpha}{Numeric reflection coefficient.  Must be positive.  Defaults to 1.0.}
  \item{betao, betai}{Numeric outside- and inside-contraction coefficients, respectively.  Both must be within unit interval (0,1).  Both default to 0.5.}
  \item{gamma}{Numeric expansion coefficient.  If positive, must be greater than \code{alpha}.  If non-positive, expansion transformations will not be carried out.  Defaults to 2.0.}
  \item{sigma}{Numeric shrink coefficient.  Cannot exceed 1.0.  If non-positive, shrink coefficients will not be carried out, and failed contractions will instead be followed by a simplex restart.  Defaults to 0.5.}
  \item{bignum}{Numeric value with which the fitfunction value is to be replaced if the fit is non-finite or is evaluated at infeasible parameter values.  Defaults to 1e35.}
  \item{iniSimplexType}{Character string naming the method by which to construct the initial simplex from the free-parameter start values.  Defaults to "regular".}
  \item{iniSimplexEdge}{Numeric edge-length of the initial simplex.  Defaults to 1.0.}
  \item{iniSimplexMat}{Optional numeric matrix providing the vertices of the initial simplex.  The matrix must have as many columns as there are free parameters in the \link{MxModel}.  The matrix's number of rows must be no less than the number of free parameters minus the number of degrees-of-freedom gained from equality \link[mxConstraint]{constraints}, if any.  If a non-\code{NA} value is provided, argument \code{iniSimplexEdge} is ignored, and argument \code{iniSimplexType} is only used in the case of a restart.}
  \item{greedyMinimize}{Logical; should the optimizer use "greedy minimization?"  Defaults to \code{FALSE}.  See below for details.}
  \item{altContraction}{Logical; should the optimizer use an "alternate contraction" transformation?  Defaults to \code{FALSE}.  See below for details.}
  \item{degenLimit}{Numeric "degeneracy limit;" defaults to 0.  If positive, the simplex will be restarted if the measure of the angle between any two of its edges is within 0 or pi by less than \code{degenLimit}.}
  \item{stagnCtrl}{"Stagnation control;" integer vector of length 2; defaults to \code{c(-1L,-1L)}.  See below for details.}
  %%\item{validationRestart}
  \item{xTolProx}{Numeric "domain-convergence" criterion; defaults to 1e-4.  See below for details.}
  \item{fTolProx}{Numeric "range-convergence" criterion; defaults to 1e-4.  See below for details.}
  %%\item{doPseudoHessian}
  \item{ineqConstraintMthd}{"Inequality constraint method;" character string. Defaults to "soft", which is the only such method currently implemented.}
  \item{eqConstraintMthd}{"Equality constraint method;" character string. Defaults to "soft", which is the only such method currently implemented.}
}
\details{
An iteration of the original (1965) Nelder-Mead algorithm proceeds as follows (only displayed in the PDF format of this documentation):
\if{latex}{
	%%TODO
}


%%other details TODO

}
\value{
Returns an object of class '\link{MxComputeNelderMead}'.
}

\references{
Gill, P. E., Murray, W., & Wright, M. H.  (1982).  \emph{Practical Optimization}.  Bingley, UK: Emerald Group Publishing Ltd.

Kelley, C. T.  (1999).  Detection and remediation of stagnation in the Nelder-Mead algorithm using a sufficient decrease condition.  \emph{SIAM Journal of Optimization 10}(1), p. 43-55.

Nelder, J. A., & Mead, R.  (1965) .  A simplex method for function minimization.  \emph{The Computer Journal, 7}, p. 308-313.

Singer, S., & Nelder, J.  (2009).  Nelder-Mead algorithm.  \emph{Scholarpedia, 4}(7):2928., revision #91557.  http://www.scholarpedia.org/article/Nelder-Mead_algorithm .

Singer, S., & Singer, S.  (2004).  Efficient implementation of the Nelder-Mead search algorithm.  \emph{Applied Numerical Analysis & Computational Mathematics Journal, 1}(2), p. 524-534.  doi: 10.1002/anac.200410015
}

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}

\examples{
foo <- mxComputeNelderMead()
str(foo)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
