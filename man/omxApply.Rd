\name{omxApply}
\alias{omxApply}

\title{On-Demand Parallel Apply}

\description{
   If the snowfall library is loaded, then this function calls \code{\link[snowfall]{sfApply}}.  Otherwise it invokes \code{\link{apply}}.
}

\usage{
omxApply(x, margin, fun, ...)
}

\arguments{
   \item{x}{a vector (atomic or list) or an expressions vector. Other objects (including classed objects) will be coerced by \code{\link{as.list}}.}
   \item{margin}{a vector giving the subscripts which the function will be applied over.}
   \item{fun}{the function to be applied to each element of \code{x}.}
   \item{...}{optional arguments to \code{fun}.}
}

\seealso{
\code{\link{omxLapply}}, \code{\link{omxSapply}} 
}

\examples{

x <- cbind(x1 = 3, x2 = c(4:1, 2:5))
dimnames(x)[[1]] <- letters[1:8]
omxApply(x, 2, mean, trim = .2)

}
