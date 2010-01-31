\name{omxSapply}
\alias{omxSapply}

\title{On-Demand Parallel Sapply}

\description{
   If the snowfall library is loaded, then this function calls \code{\link[snowfall]{sfSapply}}.  Otherwise it invokes \code{\link{sapply}}.
}

\usage{
omxSapply(x, fun, ..., simplify = TRUE, USE.NAMES = TRUE)
}

\arguments{
   \item{x}{a vector (atomic or list) or an expressions vector. Other objects (including classed objects) will be coerced by \code{\link{as.list}}.}
   \item{fun}{the function to be applied to each element of \code{x}.}
   \item{...}{optional arguments to \code{fun}.}
   \item{simplify}{logical; should the result be simplified to a vector or matrix if possible?}
   \item{USE.NAMES}{logical; if \code{TRUE} and if \code{x} is a character, use \code{x} as \code{\link{names}} for the result unless it had names already.}
}

\seealso{
\code{\link{omxApply}}, \code{\link{omxLapply}} 
}

\examples{

x <- list(a = 1:10, beta = exp(-3:3), logic = c(TRUE,FALSE,FALSE,TRUE))
# compute the list mean for each list element
omxSapply(x, quantile)

}
