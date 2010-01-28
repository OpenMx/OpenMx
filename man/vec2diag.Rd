\name{vec2diag}
\alias{vec2diag}

\title{Create Diagonal Matrix From Vector}

\description{
   Given an input row or column vector, \code{vec2diag} returns a diagonal matrix with the input argument along the diagonal.
}

\usage{
vec2diag(x)
}

\arguments{
   \item{x}{a row or column vector.}
}

\details{
   Similar to the function \code{\link{diag}}, except that the input argument is always
   treated as a vector of elements to place along the diagonal.
}

\seealso{
\code{\link{diag2vec}} 
}

\examples{

vec2diag(matrix(1:4, 1, 4))
vec2diag(matrix(1:4, 4, 1))


}
