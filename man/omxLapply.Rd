\name{omxLapply}
\alias{omxLapply}

\title{On-Demand Parallel Lapply}

\description{
   If the snowfall library is loaded, then this function calls \code{\link[snowfall]{sfLapply}}.  Otherwise it invokes \code{\link{lapply}}.
}

\usage{
omxLapply(x, fun, ...)
}

\arguments{
   \item{x}{a vector (atomic or list) or an expressions vector. Other objects (including classed objects) will be coerced by \code{\link{as.list}}.}
   \item{fun}{the function to be applied to each element of \code{x}.}
   \item{...}{optional arguments to \code{fun}.}
}

\seealso{
\code{\link{omxApply}}, \code{\link{omxSapply}} 
}

\examples{

x <- list(a = 1:10, beta = exp(-3:3), logic = c(TRUE,FALSE,FALSE,TRUE))
# compute the list mean for each list element
omxLapply(x,mean)

}
