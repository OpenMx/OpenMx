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
The mxAlgebra function is used to create algebraic expressions that operate on one or more \link{MxMatrix} objects. To evaluate an \link{MxAlgebra} object, it must be placed in an \link{MxModel} object, along with all referenced \code{MxMatrix} objects and the \code{mxAlgebraObjective} function. The \code{mxAlgebraObjective} function must reference the \code{MxAlgebra} object to be evaluated by name.

The following operations are supported in mxAlgebra:
\itemize{
\item{\code{solve()}}{Inversion}
\item{\code{t()}}{Transposition}
\item{\code{+}}{Addition}
\item{\code{-}}{Subtraction}
\item{\code{\%*\%}}{Matrix Multiplication}
\item{\code{*}}{Element or dot product}
\item{\code{/}}{Element division}
\item{\code{\%x\%}}{Kronecker product}
\item{\code{\%&\%}}{Quadratic product}
\item{\code{cbind()}}{Horizontal adhesion}
\item{\code{rbind()}}{Vertical adhesion}
}

The following functions are supported in mxAlgebra:
\itemize{
\item{\code{sum(diag())}}{Trace}
\item{\code{det()}}{Determinant}
\item{\code{sum()}}{Sum}
\item{\code{max()}}{Maximum}
\item{\code{min()}}{Min}
\item{\code{abs()}}{Absolute value}
\item{\code{sin()}}{Sine}
\item{\code{sinh()}}{Hyperbolic sine}
\item{\code{cos()}}{Cosine}
\item{\code{cosh()}}{Hyperbolic cosine}
\item{\code{tan()}}{Tangent}
\item{\code{tanh()}}{Hyperbolic tangent}
\item{\code{exp()}}{Exponent}
\item{\code{log()}}{Natural Logarithm}
\item{\code{sqrt()}}{Square root}
\item{\code{diag()}}{Diagonal to vector or vector to diagonal, depending on input}
\item{\code{c(t()}}{Matrix to vector (by row; returns vector)}
\item{\code{c()}}{Matrix to vector (by column)}
\item{\code{vech()}}{Lower triangular matrix to vector (requires MCMCpack package)}
\item{\code{cov2cor()}}{Standardize matrix}
\item{\code{Re(eigen(A, only.values=TRUE)$values)}}{Real eigenvalues (for matrix A)}
\item{\code{Re(eigen(A, )$values)}}{Real eigenvectors (for matrix A)}
\item{\code{Im(eigen(A, only.values=TRUE)$values)}}{Imanginary eigenvalues (for matrix A)}
\item{\code{Im(eigen(A, )$values)}}{Imanginary eigenvectors (for matrix A)}
\item{\code{colMeans()}}{Mean of columns}
\item{\code{var()}}{Covariance of columns}
}}

\value{
	Returns a new \link{MxAlgebra} object.
}

\references{
Algebra operations and discussion of supported functions may be found in Chapter 4 of the Mx manual, which can be found at http://www.vcu.edu/mx/documentation.html.
}

\seealso{
\link{MxAlgebra} for the S4 class created by mxAlgebra. \link{MxMatrix} and \link{mxMatrix} for objects which may be entered in the 'expression' argument and the function that creates them. More information about the OpenMx package may be found \link[=OpenMx]{here}. 
}

\examples{

# the algebra B evaluates to the matrix A
A <- mxMatrix("Full", nrow = 3, ncol = 3, name = "A")
B <- mxAlgebra(A, name = "B")

# the algebra C is A + B
C <- mxAlgebra(A + B, name = "C")

\dontrun{ % not all functions have been implemented yet
# the algebra D is sin(C)
D <- mxAlgebra(sin(C), name = "D")

# Evaluate the mxAlgebra object 'D'
modelD <- mxModel(model="AlgebraExample", A, B, C, D, mxAlgebraObjective("D"))
modelD <- mxRun(modelD)
}

}
