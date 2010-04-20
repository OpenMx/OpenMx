\name{eigenvec}
\alias{eigenval}
\alias{eigenvec}
\alias{ieigenval}
\alias{ieigenvec}

\title{Eigenvector/Eigenvalue Decomposition}

\description{
   \code{eigenval} computes the real parts of the eigenvalues of a square matrix.
   \code{eigenval} computes the real parts of the eigenvectors of a square matrix.
   \code{ieigenval} computes the imaginary parts of the eigenvalues of a square matrix.
   \code{ieigenvec} computes the imaginary parts of the eigenvectors of a square matrix.
   \code{eigenval} and \code{ieigenval} return nx1 matrices containing the real or imaginary parts of the eigenvalues, sorted in decreasing order of the modulus of the complex eigenvalue.  For eigenvalues without an imaginary part, this is equivalent to sorting in decreasing order of the absolute value of the eigenvalue. (See \code{\link{Mod}} for more info.)
   \code{eigenvec} and \code{ieigenvec} return nxn matrices, where each column corresponds to an eigenvector.  These are sorted in decreasing order of the modulus of their associated complex eigenvalue.
}

\usage{
eigenval(x)
eigenvec(x)
ieigenval(x)
ieigenvec(x)
}

\arguments{
   \item{x}{the square matrix whose eigenvalues/vectors are to be calculated.}
}

\details{

    Eigenvectors returned by \code{eigenvec} and \code{ieigenvec} are normalized to unit length.

}

\seealso{
\code{\link{eigen}}
}

\examples{

A <- mxMatrix(values = runif(25), nrow = 5, ncol = 5, name = 'A')
G <- mxMatrix(values = c(0, -1, 1, -1), nrow=2, ncol=2, name='G')

model <- mxModel(A, G, name = 'model')

mxEval(eigenvec(A), model)
mxEval(eigenvec(G), model)
mxEval(eigenval(A), model)
mxEval(eigenval(G), model)
mxEval(ieigenvec(A), model)
mxEval(ieigenvec(G), model)
mxEval(ieigenval(A), model)
mxEval(ieigenval(G), model)



}
