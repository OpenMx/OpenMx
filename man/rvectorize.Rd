\name{rvectorize}
\alias{rvectorize}

\title{Vectorize By Row}

\description{
   This function returns the vectorization of an input matrix in a row by row traversal of the matrix.  The output is returned as a column vector.
}

\usage{
rvectorize(x)
}

\arguments{
   \item{x}{an input matrix.}
}

\seealso{
\code{\link{cvectorize}}, \code{\link{vech}}, \code{\link{vechs}} 
}

\examples{

rvectorize(matrix(1:9, 3, 3))
rvectorize(matrix(1:12, 3, 4))


}
