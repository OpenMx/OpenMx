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
	iniSimplexEdge=1, iniSimplexMat=NULL, greedyMinimize=FALSE, 
	altContraction=FALSE, degenLimit=0, stagnCtrl=c(-1L,-1L),
	validationRestart=TRUE,
	xTolProx=1e-8, fTolProx=1e-8,
	doPseudoHessian=TRUE,
	ineqConstraintMthd=c("soft","eqMthd"), 
	eqConstraintMthd=c("GDsearch","soft","backtrack","l1p"),
	backtrackCtrl=c(0.5,5),
	centerIniSimplex=FALSE)
}
\arguments{
  \item{freeSet}{Character-string names of \link[=MxMatrix]{MxMatrices} containing free parameters.}
  \item{fitfunction}{Character-string name of the fitfunction; defaults to 'fitfunction'.}
  \item{verbose}{Integer level of reporting printed to terminal at \link[=mxRun]{runtime}; defaults to 0.}
  \item{nudgeZeroStarts}{Should free parameters with start values of zero be "nudged" to 0.1 at \link[=mxRun]{runtime}?  Defaults to the current global value of \link{mxOption} "Nudge zero starts".  May be a logical value, or one of character strings "Yes" or "No".}
  \item{maxIter}{Integer maximum number of iterations.  Value of \code{NULL} is accepted, in which case the value used at \link[=mxRun]{runtime} will be 10 times the number of iterations specified by the effective value of \link{mxOption} "Major iterations".}
  \item{...}{Not used.  Forces remaining arguments to be specified by name.}
  \item{alpha}{Numeric reflection coefficient.  Must be positive.  Defaults to 1.0.}
  \item{betao, betai}{Numeric outside- and inside-contraction coefficients, respectively.  Both must be within unit interval (0,1).  Both default to 0.5.}
  \item{gamma}{Numeric expansion coefficient.  If positive, must be greater than \code{alpha}.  If non-positive, expansion transformations will not be carried out.  Defaults to 2.0.}
  \item{sigma}{Numeric shrink coefficient.  Cannot exceed 1.0.  If non-positive, shrink transformations will not be carried out, and failed contractions will instead be followed by a simplex restart.  Defaults to 0.5.}
  \item{bignum}{Numeric value with which the fitfunction value is to be replaced if the fit is non-finite or is evaluated at infeasible parameter values.  Defaults to 1e35.}
  \item{iniSimplexType}{Character string naming the method by which to construct the initial simplex from the free-parameter start values.  Defaults to "regular".}
  \item{iniSimplexEdge}{Numeric edge-length of the initial simplex.  Defaults to 1.0.}
  \item{iniSimplexMat}{Optional numeric matrix providing the vertices of the initial simplex.  The matrix must have as many columns as there are free parameters in the \link{MxModel}.  The matrix's number of rows must be no less than the number of free parameters minus the number of degrees-of-freedom gained from equality \link[=mxConstraint]{MxConstraints}, if any.  If a non-\code{NULL} value is provided, argument \code{iniSimplexEdge} is ignored, and argument \code{iniSimplexType} is only used in the case of a restart.}
  \item{greedyMinimize}{Logical; should the optimizer use "greedy minimization?"  Defaults to \code{FALSE}.  See below for details.}
  \item{altContraction}{Logical; should the optimizer use an "alternate contraction" transformation?  Defaults to \code{FALSE}.  See below for details.}
  \item{degenLimit}{Numeric "degeneracy limit;" defaults to 0.  If positive, the simplex will be restarted if the measure of the angle between any two of its edges is within 0 or pi by less than \code{degenLimit}.}
  \item{stagnCtrl}{"Stagnation control;" integer vector of length 2; defaults to \code{c(-1L,-1L)}.  See below for details.}
  \item{validationRestart}{Logical; defaults to TRUE.}
  \item{xTolProx}{Numeric "domain-convergence" criterion; defaults to 1e-8.  See below for details.}
  \item{fTolProx}{Numeric "range-convergence" criterion; defaults to 1e-8.  See below for details.}
  \item{doPseudoHessian}{Logical; defaults to \code{TRUE}.}
  \item{ineqConstraintMthd}{"Inequality constraint method;" character string. Defaults to "soft".}
  \item{eqConstraintMthd}{"Equality constraint method;" character string. Defaults to "GDsearch".}
  \item{backtrackCtrl}{Numeric vector of length two.  See below for details.}
  \item{centerIniSimplex}{Logical. If \code{FALSE} (default), the MxModel's start values are used as the "first" vertex of the initial simplex.  If \code{TRUE}, the initial simplex is re-centered so that the MxModel's start values are its eucentroid.  However, if \code{iniSimplexMat} is non-\code{NULL} or if \code{iniSimplexType="smartRight"}, a value of \code{TRUE} is treated as \code{FALSE}.}
}
\details{
%TODO: embellish with latex

The state of a Nelder-Mead optimization problem is represented by a simplex (polytope) of \eqn{n+1} vertices in the space of the free parameters, where \eqn{n} is the number of free parameters minus the number of degrees-of-freedom gained from equality \link[=mxConstraint]{MxConstraints}.  An iteration of the algorithm first sorts the \eqn{n+1} vertices by their corresponding fitfunction values (i.e., the values of the fitfunction when evaluated at each vertex), in ascending order (i.e., from "best" fit to "worst" fit).  Then, the "subcentroid," which is the centroid of the "best" \eqn{n} vertices, is calculated.  Then, the algorithm attempts to improve upon the worst fit by transforming the simplex; see \href{http://www.scholarpedia.org/article/Nelder-Mead_algorithm}{Singer & Nelder (2009)} for details.

Argument \code{iniSimplexType} dictates how the initial simplex will be constructed from the start values if argument \code{iniSimplexMat} is \code{NULL}, and how the simplex will be re-initialized in the case of a restart.  In all four cases, the vector of start values constitutes the "starting vertex" of the initial simplex.  If \code{iniSimplexType="regular"}, the initial simplex is merely a regular simplex with edge length equal to \code{iniSimplexEdge}.  A \code{"right"} simplex is constructed by incrementing each free parameter by \code{iniSimplexEdge} from its starting value; thus, all the edges that intersect at the starting vertex do so at right angles.  A \code{"smartRight"} simplex is constructed similarly, except that each free parameter is both incremented \emph{and} decremented by \code{iniSimplexEdge}, and of those two points the one with the smaller fitfunction value is retained as a vertex.  A \code{"random"} simplex is constructed by randomly perturbing the start values, in a manner similar to the default for \code{\link{mxTryHard}()}, to generate the coordinates of the other vertices.  The user is advised that bounds on the free parameters may keep the initial simplex from having the requested regularity or edge-length, and that \code{iniSimplexType} is at best a \emph{suggestion} in the presence of equality \link[=mxConstraint]{MxConstraints}.

Note that if argument \code{iniSimplexMat} has nonzero length, the actual start values of the MxModel's free parameters are not used as a vertex of the initial simplex (unless one of the rows of \code{iniSimplexMat} happens to contain those start values).

If the simplex is restarted, a new simplex is constructed per argument \code{iniSimplexType}, with edge length equal to the distance between the current best and second-best vertices, and with the current best vertex used as the "first" vertex.

If \code{greedyMinimize=FALSE}, "greedy expansion" (Singer & Singer, 2004) is used: if the expansion point and reflection point both have smaller fitfunction values than the best vertex, the expansion point is accepted.  If \code{greedyMinimize=TRUE}, "greedy minimization" (Singer & Singer, 2004) is used: if the expansion point and the reflection point both have smaller fitfunction values than the best vertex, the better of the two new points is accepted.

If argument \code{altContraction=TRUE}, the "modified contraction step" of Gill et al. (1982, Chapter 4) is used, and the candidate point is contracted toward the best vertex instead of toward the subcentroid.

If positive, the first element of argument \code{stagnCtrl} sets a threshold for the number of successive iterations in which the best vertex of the simplex does not change, after which the algorithm is said to be "stagnant" (in a sense similar to that of Kelley, 1999).  To attempt to remedy the stagnation, the simplex is restarted.  If positive, the second element of argument \code{stagnCtrl} sets threshold for the number of restarts conducted, beyond which stagnation no longer triggers a restart.  The rationale for the second element is that the best vertex may not change for many iterations when the optimizer is close to convergence, under which circumstances restarting would be counterproductive, and in any event would require additional fitfunction evaluations.

If argument \code{validationRestart=TRUE}, then when the optimizer has successfully converged, it will restart the simplex and attempt to improve upon the tentative solution it already found.  This validation restart (Gill et al., 1982, Chapter 4) always re-initializes the simplex as a regular simplex, centered on the best vertex of the tentative solution, with edge-length equal to the distance between the best and worst vertices of the tentative solution.  Optimization proceeds until convergence to a solution with a better fit value, or \eqn{2n} iterations have elapsed.

The Nelder-Mead optimizer is considered to have successfully converged if (1) the largest \emph{l}-infinity norm of the vector-differences between the best vertex and the other vertices is less than argument \code{xTolProx}, or (2) if the largest absolute difference in fit value between the best vertex and the other vertices is less than \code{fTolProx}.

If argument \code{doPseudoHessian=TRUE}, there are no equality \link[=mxConstraint]{MxConstraints}, and the "l1p" method (see below) is not in use for inequality \link[=mxConstraint]{MxConstraints}, then OpenMx will attempt to calculate the "pseudo-Hessian" or "curvature" matrix as described in the appendix to Nelder & Mead (1965).  If successful, this matrix will be stored in the 'output' slot of the post-\link[=mxRun]{run} MxComputeNelderMead object.  Although crude, its inverse can be used as an estimate of the repeated-sampling covariance matrix of the free parameters when the usual finite-differences Hessian is unreliable.

OpenMx's implementation of Nelder-Mead can handle nonlinear inequality \link[=mxConstraint]{MxConstraints} reasonably well.  Its default method for doing so, with argument \code{ineqConstraintMthd="soft"}, imposes a "soft" feasibility constraint by assigning a fitfunction value of \code{bignum} to points that violate the constraints by more than \link{mxOption} 'Feasibility tolerance'.  Alternately, with argument \code{ineqConstraintMthd="eqMthd"}, inequality \link[=mxConstraint]{MxConstraints} can be handled by the same method provided to argument \code{eqConstraintMthd}, whether or not equality \link[=mxConstraint]{MxConstraints} are present.

OpenMx's implementation of Nelder-Mead respects equality \link[=mxConstraint]{MxConstraints}, but does not handle them especially well.  Its effectiveness at handling equalities may be improved by providing a matrix to argument \code{iniSimplexMat} that ensures \emph{all} of the initial vertices are feasible.  Users are warned that this Nelder-Mead implementation will not work correctly with MxModels containing redundant equality MxConstraints, and presently has no way of detecting whether any are present.  If argument \code{eqConstraintMthd="GDsearch"} (the default), then whenever Nelder-Mead evaluates the fitfunction at an infeasible point, it initiates a subsidiary optimization that uses \link[=mxComputeGradientDescent]{SLSQP} to find the nearest (in squared Euclidean distance) feasible point, and replaces that feasible point for the infeasible one.  The user should note that the function evaluations that occur during this subsidiary optimization are counted toward the total number of fitfunction evaluations during the call to \code{\link{mxRun}()}.  The effectiveness of the 'GDsearch' method is often improved by setting \link{mxOption} 'Feasibility tolerance' to a stricter (smaller) value than the on-load default.  The method specified by \code{eqConstraintMthd="soft"} is described in the preceding paragraph.  If argument \code{eqConstraintMthd="backtrack"}, then the optimizer attempts to backtrack from an infeasible point to a feasible point in a manner similar to that of Ghiasi et al. (2008), except that it used with \emph{all} new points, and not just those encountered via reflection, expansion and contraction.  In this case, the displacement from the prior point to the candidate point is reduced by the proportion provided as the first element of argument \code{backtrackCtrl}, and thus a new candidate point is considered.  This process is repeated until feasibility of the candidate point is restored, or the number of attempts exceeds the second element of argument \code{backtrackCtrl}.   If argument \code{eqConstraintMthd="l1p"}, Nelder-Mead is used as part of an \eqn{l_1}{l1}-penalty algorithm.  When using "l1p", the simplex gradient (Kelley, 1999) and "pseudo-Hessian" are never calculated. 
}
\value{
Returns an object of class 'MxComputeNelderMead'.
}

\references{
Ghiasi, H., Pasini, D., & Lessard, L.  (2008).  Constrained globalized Nelder-Mead method for simultaneous structural and manufacturing optimization of a composite bracket.  \emph{Journal of Composite Materials, 42}(7), p. 717-736.  doi: 10.1177/0021998307088592

Gill, P. E., Murray, W., & Wright, M. H.  (1982).  \emph{Practical Optimization}.  Bingley, UK: Emerald Group Publishing Ltd.

Kelley, C. T.  (1999).  Detection and remediation of stagnation in the Nelder-Mead algorithm using a sufficient decrease condition.  \emph{SIAM Journal of Optimization 10}(1), p. 43-55.

Nelder, J. A., & Mead, R.  (1965) .  A simplex method for function minimization.  \emph{The Computer Journal, 7}, p. 308-313.

Singer, S., & Nelder, J.  (2009).  Nelder-Mead algorithm.  \emph{Scholarpedia, 4}(7):2928., revision #91557.  http://www.scholarpedia.org/article/Nelder-Mead_algorithm .

Singer, S., & Singer, S.  (2004).  Efficient implementation of the Nelder-Mead search algorithm.  \emph{Applied Numerical Analysis & Computational Mathematics Journal, 1}(2), p. 524-534.  doi: 10.1002/anac.200410015
}

\examples{
foo <- mxComputeNelderMead()
str(foo)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
