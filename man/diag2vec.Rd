\name{diag2vec}
\alias{diag2vec}

\title{Extract Diagonal of a Matrix}

\description{
   Given an input matrix, \code{diag2vec} returns a column vector of the elements along the diagonal.
}

\usage{
diag2vec(x)
}

\arguments{
   \item{x}{an input matrix.}
}

\details{
   Similar to the function \code{\link{diag}}, except that the input argument is always
   treated as a matrix.
}

\seealso{
\code{\link{vec2diag}} 
}

\examples{

diag2vec(matrix(1:9, 3, 3))
diag2vec(matrix(1:12, 3, 4))


}
