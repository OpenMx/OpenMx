\name{mxAlgebra}
\alias{mxAlgebra}

\title{Create MxAlgebra Object}

\description{
   This function creates a new \link{MxAlgebra} object.
}

\usage{
mxAlgebra(expression, name = NA, dimnames = NA)
}

\arguments{
   \item{expression}{An R expression of matrix operators and matrix functions.}
   \item{name}{An optional character string indicating the name of the object.}
   \item{dimnames}{list. The dimnames attribute for the algebra: a list of length 2 giving the row and column names respectively. An empty list is treated as NULL, and a list of length one as row names. The list can be named, and the list names will be used as names for the dimensions.}
}

\details{
The mxAlgebra function is used to create algebraic expressions that operate on one or more \link{MxMatrix} objects. To evaluate an \link{MxAlgebra} object, it must be placed in an \link{MxModel} object, along with all referenced \code{MxMatrix} objects and the \code{mxAlgebraObjective} function. The \code{mxAlgebraObjective} function must reference the \code{MxAlgebra} object to be evaluated by name.

The following operations are supported in mxAlgebra:
\describe{
\item{\code{solve()}}{Inversion}
\item{\code{t()}}{Transposition}
\item{\code{^}}{Elementwise powering}
\item{\code{\%^\%}}{Kronecker powering}
\item{\code{+}}{Addition}
\item{\code{-}}{Subtraction}
\item{\code{\%*\%}}{Matrix Multiplication}
\item{\code{*}}{Element or dot product}
\item{\code{/}}{Element division}
\item{\code{\%x\%}}{Kronecker product}
\item{\code{\%&\%}}{Quadratic product}
}

The following functions are supported in mxAlgebra:
\describe{
\item{\code{cbind}}{Horizontal adhesion}
\item{\code{rbind}}{Vertical adhesion}
\item{\code{det}}{Determinant}
\item{\code{tr}}{Trace}
\item{\code{sum}}{Sum}
\item{\code{prod}}{Product}
\item{\code{max}}{Maximum}
\item{\code{min}}{Min}
\item{\code{abs}}{Absolute value}
\item{\code{sin}}{Sine}
\item{\code{sinh}}{Hyperbolic sine}
\item{\code{cos}}{Cosine}
\item{\code{cosh}}{Hyperbolic cosine}
\item{\code{tan}}{Tangent}
\item{\code{tanh}}{Hyperbolic tangent}
\item{\code{exp}}{Exponent}
\item{\code{log}}{Natural Logarithm}
\item{\code{sqrt}}{Square root}
\item{\code{\link{rvectorize}}}{Vectorize by row}
\item{\code{\link{cvectorize}}}{Vectorize by column}
\item{\code{\link{vech}}}{Half-vectorization}
\item{\code{\link{vechs}}}{Strict half-vectorization}
\item{\code{\link{vec2diag}}}{Create a diagonal matrix}
\item{\code{\link{diag2vec}}}{Extract diagonal from matrix}
}}

\value{
	Returns a new \link{MxAlgebra} object.
}

\references{
The OpenMx User's guide can be found at http://openmx.psyc.virginia.edu/documentation.
}

\seealso{
\link{MxAlgebra} for the S4 class created by mxAlgebra. \link{MxMatrix} and \link{mxMatrix} for objects which may be entered in the 'expression' argument and the function that creates them. More information about the OpenMx package may be found \link[=OpenMx]{here}. 
}

\examples{

A <- mxMatrix("Full", nrow = 3, ncol = 3, values=2, name = "A")

# Simple example: algebra B simply evaluates to the matrix A
B <- mxAlgebra(A, name = "B")

# Compute A + B
C <- mxAlgebra(A + B, name = "C")

# Compute sin(C)
D <- mxAlgebra(sin(C), name = "D")

# Make a model and evaluate the mxAlgebra object 'D'
A <- mxMatrix("Full", nrow = 3, ncol = 3, values=2, name = "A")
model <- mxModel("AlgebraExample", A, B, C, D )
fit   <- mxRun(model)
mxEval(D, fit)

}
