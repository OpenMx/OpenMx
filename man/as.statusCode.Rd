\name{as.statusCode}
\alias{as.statusCode}

\title{Convert a numeric or character vector into an optimizer status code factor}

\usage{
as.statusCode(code)
}

\arguments{
   \item{code}{a character or numeric vector of optimizer status code}
 }

\description{

   \itemize{
\item 0,\sQuote{OK}: Optimization succeeded
   \item 1,\sQuote{OK/green}: Optimization succeeded, but the sequence
   of iterates has not yet converged (Mx status GREEN). This condition is only detected by
   NPSOL.
   \item 2,\sQuote{infeasible linear constraint}:
   The linear constraints and bounds could not be satisfied. The problem has no feasible solution.
   \item 3,\sQuote{infeasible non-linear constraint}:
   The nonlinear constraints and bounds could not be satisfied. The problem may have no feasible solution.
   \item 4,\sQuote{iteration limit}:
   Optimization was stopped prematurely because the iteration limit was
   reached (Mx status BLUE). You might want to rerun: \code{m1 = mxRun(m1)} or increase the
   iteration limit (see \code{\link{mxOption}}).
   \item 5,\sQuote{not convex}:
   The Hessian at the solution does not appear to be convex (Mx status
   RED). There may be more than one solution to the model. See
   \code{\link{mxCheckIdentification}}.
   \item 6,\sQuote{nonzero gradient}:
   The model does not satisfy the first-order optimality conditions to
   the required accuracy, and no improved point for the merit function
   could be found during the final linesearch (Mx status RED).
   To search nearby, see \code{\link{mxTryHard}}.
   \item 7,\sQuote{bad deriv}:
   You have provided analytic derivatives. However,
   your provided derivatives differ too much from numerically
   approximated derivatives. Double check your math.
   \item 9,\sQuote{internal error}: An input parameter was invalid. The
   most likely cause is a bug in the code. Please report occurrences to
   the OpenMx developers.
   \item 10,\sQuote{infeasible start}:
   Starting values were infeasible.
   Modify the start values for one or more parameters. For instance,
   set means to their measured value, or set variances and covariances to
   plausible values. See \code{\link{mxAutoStart}} and \code{\link{mxTryHard}}.
   }
}

\seealso{
\code{\link{mxBootstrap}}
\code{\link{summary.MxModel}}
}
