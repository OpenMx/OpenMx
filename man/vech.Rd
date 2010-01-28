\name{vech}
\alias{vech}

\title{Half-vectorization}

\description{
   This function returns the half-vectorization of an input matrix as a column vector.
}

\usage{
vech(x)
}

\arguments{
   \item{x}{an input matrix.}
}

\details{
   The half-vectorization of an input matrix consists of the elements in the lower triangle of the matrix, including the elements along the diagonal of the matrix, as a column vector.  The column vector is created by traversing the matrix in column-major order.
}

\seealso{
\code{\link{vech}} 
}

\examples{

vech(matrix(1:9, 3, 3))
vech(matrix(1:12, 3, 4))


}
