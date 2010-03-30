\name{vechs}
\alias{vechs}

\title{Strict Half-vectorization}

\description{
   This function returns the strict half-vectorization of an input matrix as a column vector.
}

\usage{
vechs(x)
}

\arguments{
   \item{x}{an input matrix.}
}

\details{
   The half-vectorization of an input matrix consists of the elements in the lower triangle of the matrix, excluding the elements along the diagonal of the matrix, as a column vector.  The column vector is created by traversing the matrix in column-major order.
}

\seealso{
\code{\link{vech}}, \code{\link{rvectorize}}, \code{\link{cvectorize}}
}

\examples{

vechs(matrix(1:9, 3, 3))
vechs(matrix(1:12, 3, 4))


}
