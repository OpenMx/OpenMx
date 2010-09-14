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
   treated as a matrix (i.e., it doesn't have diag()'s functions of returning an Identity matrix from an nrow specification, nor to return 
   a matrix wrapped around a diagonal if provided with a vector). To get vector2matrix functionality, call \code{\link{vec2diag}}.
}

\seealso{
\code{\link{vec2diag}} 
}

\examples{

diag2vec(matrix(1:9, nrow=3))
#      [,1]
# [1,]    1
# [2,]    5
# [3,]    9

diag2vec(matrix(1:12, nrow=3, ncol=4))
#      [,1]
# [1,]    1
# [2,]    5
# [3,]    9

}
