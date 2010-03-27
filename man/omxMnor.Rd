\name{omxMnor}
\alias{omxMnor}

\title{Multivariate Normal Integration}

\description{
   Given a covariance matrix, a means vector, and vectors of lower and upper bounds, returns the multivariate normal integral across the space between bounds.
}

\usage{
omxMnor(covariance, means, lbound, ubound)
}

\arguments{
   \item{covariance}{the covariance matrix describing the multivariate normal distribution.}
   \item{means}{a row vector containing means of the variables of the underlying distribution.}
   \item{lbound}{a row vector containing the lower bounds of the integration in each variable.}
   \item{ubound}{a row vector containing the upper bounds of the integration in each variable.}
}

\details{
   The order of columns in the \sQuote{means}, \sQuote{lbound}, and \sQuote{ubound} vectors are assumed to be the same as that of the covariance matrix.  That is, means[i] is considered to be the mean of the variable whose variance is in covariance[i,i].  That variable will be integrated from lbound[i] to ubound[i] as part of the integration.

   The value of ubound[i] or lbound[i] may be set to Inf or -Inf if a boundary at positive or negative infinity is desired.
   
   For all i, ubound[i] must be strictly greater than lbound[i].
}

\examples{

data(myFAData)

covariance <- cov(myFAData[,1:3])
means <- mean(myFAData[,1:3])
lbound <- c(-Inf, 0,   1)    # Integrate from -Infinity to 0 on first variable
ubound <- c(0,    Inf, 2.5)  # From 0 to +Infinity on second, and from 1 to 2.5 on third
omxMnor(covariance, means, lbound, ubound)

}
