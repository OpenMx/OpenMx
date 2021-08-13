\name{omxDefaultComputePlan}
\alias{omxDefaultComputePlan}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Construct default compute plan
}
\description{
This function generates a default compute plan, where "default compute plan" refers to an object of class \link[=mxComputeSequence]{MxComputeSequence} which is appropriate for use in a wide variety of cases.  The exact specification of the plan will depend upon the arguments provided to \code{omxDefaultComputePlan()}.
}
\usage{
omxDefaultComputePlan(modelName=NULL, intervals=FALSE,
  useOptimizer=TRUE, optionList=options()$mxOption)
}

\arguments{
  \item{modelName}{Optional (defaults to \code{NULL}) character string, providing the name of the \link{MxModel} the fitfunction of which is to be evaluated, and usually, optimized.}
  \item{intervals}{Logical; will \link[=mxCI]{confidence intervals} be computed? Defaults to \code{FALSE}.}
  \item{useOptimizer}{Logical; will a fitfunction be minimized?  Defaults to \code{TRUE}.}
  \item{optionList}{List of \link{mxOption}s.  Defaults to the current list of global \link{mxOption}s.}
}

\details{
At minimum, argument \code{optionList} must include \dQuote{Gradient algorithm}, \dQuote{Gradient iterations}, \dQuote{Gradient step size}, \dQuote{Calculate Hessian}, and \dQuote{Standard Errors}.
}
\value{
Returns an object of class \link[=mxComputeSequence]{MxComputeSequence}.
}

\examples{
foo <- omxDefaultComputePlan(modelName="bar")
str(foo)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
