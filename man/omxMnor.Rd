\name{omxMnor}
\alias{omxMnor}

\title{Multivariate Normal Integration}

\description{
   Given a covariance matrix, a means vector, and vectors of lower and upper bounds, returns the multivariate normal integral across the space between bounds.
}

\usage{
omxMnor(cov, means, lbounds, ubounds)
}

\arguments{
   \item{cov}{the covariance matrix describing the multivariate normal distribution.}
   \item{means}{a row vector containing means of the variables of the underlying distribution.}
   \item{lbounds}{a row vector containing the lower bounds of the integration in each variable.}
   \item{ubounds}{a row vector containing the upper bounds of the integration in each variable.}
}

\details{
   The order of columns in the means, lbounds, and ubounds vector are assumed to be the same as that of the covariance matrix.  That is, means[i] is considered to be the mean of the variable whose variance is in cov[i,i].  That variable will be integrated from lbounds[i] to ubounds[i] as part of the integration.

   The value of ubound[i] or lbound[i] may be set to Inf or -Inf if a boundary at positive or negative infinity is desired.
   
   For all i, ubound[i] must be strictly greater than lbound[i].
}

\seealso{
}

\examples{

cov <- cov(myData[,1:3])
means <- means(myData[,1:3])
lbounds <- c(-Inf, 0,   1)    % Integrate from -Infinity to 0 on first variable
ubounds <- c(0,    Inf, 2.5)   % From 0 to +Infinity on second, and from 1 to 2.5 on third
omxMnor(cov, means, lbounds, ubounds)

}
