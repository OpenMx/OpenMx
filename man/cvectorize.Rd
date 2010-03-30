\name{cvectorize}
\alias{cvectorize}

\title{Vectorize By Column}

\description{
   This function returns the vectorization of an input matrix in a column by column traversal of the matrix.  The output is returned as a column vector.
}

\usage{
cvectorize(x)
}

\arguments{
   \item{x}{an input matrix.}
}

\seealso{
\code{\link{rvectorize}}, \code{\link{vech}}, \code{\link{vechs}} 
}

\examples{

cvectorize(matrix(1:9, 3, 3))
cvectorize(matrix(1:12, 3, 4))


}
