\name{mxAlgebra}
\alias{mxAlgebra}

\title{Function To Create MxAlgebra Object}

\description{
   This function creates a new \link{MxAlgebra} object.
}

\usage{
mxAlgebra(expression, name = NA)
}

\arguments{
   \item{expression}{An R expression of matrix operators and matrix functions.}
   \item{name}{An optional character string indicating the name of the object.}
}

\details{
The mxAlgebra function is used to create algebraic expressions that operate on one or more \link{MxMatrix} objects. To evaluate an \link{MxAlgebra} object, it must be placed in an \link{MxModel} object, along with all referenced \link{MxMatrix} objects and the \link{MxAlgebraObjective} function. The \link{MxAlgebraObjective} function must reference the \link{MxAlgebra} object to be evaluated by name.

The following operations are supported in mxAlgebra:
\item{solve()} Inversion
\item{t()} Transposition
\item{+} Addition
\item{-} Subtraction
\item{/%*/%} Matrix Multiplication
\item{*} Element or dot product
\item{/} Element division
\item{/%x/%} Kronecker product
\item{/%&/%} Quadratic product
\item{cbind()} Horizontal adhesion
\item{rbind()} Vertical adhesion

The following functions are supported in mxAlgebra:
\item{sum(diag())} Trace
\item{det()} Determinant
\item{sum()} Sum
\item{max()} Maximum
\item{min()} Min
\item{abs()} Absolute value
\item{sin()} Sine
\item{sinh()} Hyperbolic sine
\item{cos()} Cosine
\item{cosh()} Hyperbolic cosine
\item{tan()} Tangent
\item{tanh()} Hyperbolic tangent
\item{exp()} Exponent
\item{log()} Natural Logarithm
\item{sqrt()} Square root
\item{diag()} Diagonal to vector or vector to diagonal, depending on input
\item{c(t()} Matrix to vector (by row; returns vector)
\item{c()} Matrix to vector (by column)
\item{vech()} Lower triangular matrix to vector (requires MCMCpack package)
\item{cov2cor()} Standardize matrix
\item{Re(eigen(A, only.values=TRUE)$values)} Real eigenvalues (for matrix A)
\item{Re(eigen(A, )$values)} Real eigenvectors (for matrix A)
\item{Im(eigen(A, only.values=TRUE)$values)} Imanginary eigenvalues (for matrix A)
\item{Im(eigen(A, )$values)} Imanginary eigenvectors (for matrix A)
\item{colMeans()} Mean of columns
\item{var()} Covariance of columns
}

\value{
	Returns a new \link{MxAlgebra} object.
}

\references{
Algebra operations and discussion of supported functions may be found in Chapter 4 of the Mx manual, which can be found at http://www.vcu.edu/mx/documentation.html.
}

\examples{

# the algebra B evaluates to the matrix A
A <- mxMatrix("Full", nrow = 3, ncol = 3, name = "A")
B <- mxAlgebra(A, name = "B")

# the algebra C is A + B
C <- mxAlgebra(A + B, name = "C")

# the algebra D is sin(C)
D <- mxAlgebra(sin(C), name = "D")

# Evaluate the mxAlgebra object 'D'
modelD<-model(model="AlgebraExample", A, B, C, D, mxAlgebraObjective("D"))
modelD<-mxJobRun(modelD)
}
