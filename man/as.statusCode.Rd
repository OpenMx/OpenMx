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
   \item 1,\sQuote{OK/green}: Optimization succeeded, but the sequence of iterates has not yet converged. NPSOL was terminated because no further improvement could be made in the merit function (Mx status GREEN)
   \item 2,\sQuote{infeasible linear constraint}:
   The linear constraints and bounds could not be satisfied. The problem has no feasible solution.
   \item 3,\sQuote{infeasible non-linear constraint}:
   The nonlinear constraints and bounds could not be satisfied. The problem may have no feasible solution.
   \item 4,\sQuote{iteration limit}:
   The major iteration limit was reached (Mx status BLUE).
   \item 5,\sQuote{not convex}:
   The Hessian at the solution does not appear to be convex (Mx status RED).
   \item 6,\sQuote{nonzero gradient}:
   The model does not satisfy the first-order optimality conditions to the required accuracy, and no improved point for the merit function could be found during the final linesearch (Mx status RED)
   \item 7,\sQuote{bad deriv}:
   The function derivates returned by funcon or funobj seem incorrect.
   \item 9,\sQuote{internal error}: An input parameter was invalid
   \item 10,\sQuote{infeasible start}:
   Starting values were infeasible.
   }
}

\seealso{
\code{\link{mxBootstrap}}
\code{\link{summary.MxModel}}
}
