%
%   Copyright 2007-2015 The OpenMx Project
%
%   Licensed under the Apache License, Version 2.0 (the "License");
%   you may not use this file except in compliance with the License.
%   You may obtain a copy of the License at
% 
%        http://www.apache.org/licenses/LICENSE-2.0
% 
%   Unless required by applicable law or agreed to in writing, software
%   distributed under the License is distributed on an "AS IS" BASIS,
%   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
%   See the License for the specific language governing permissions and
%   limitations under the License.

\name{mxAlgebra}
\alias{mxAlgebra}
\alias{p2z}
\alias{\%&\%}
\alias{\%^\%}

\title{Create MxAlgebra Object}

\description{
   This function creates a new \link{MxAlgebra} object.
}

\usage{
mxAlgebra(expression, name = NA, dimnames = NA, ..., fixed = FALSE)
}

\arguments{
   \item{expression}{An R expression of OpenMx-supported matrix operators and matrix functions.}
   \item{name}{An optional character string indicating the name of the object.}
   \item{dimnames}{list. The dimnames attribute for the algebra: a list
   of length 2 giving the row and column names respectively. An empty
   list is treated as NULL, and a list of length one as row names. The
   list can be named, and the list names will be used as names for the
   dimensions.}
   \item{...}{Not used.  Forces argument \sQuote{fixed} to be specified by name.}
   \item{fixed}{If TRUE, this algebra will not be recomputed
     automatically when things it depends on change. \link{mxComputeOnce}
     can be used to force it to recompute.}
}

\details{
The mxAlgebra function is used to create algebraic expressions that operate on one or more 
\link{MxMatrix} objects. To evaluate an \link{MxAlgebra} object, 
it must be placed in an \link{MxModel} object, along with all referenced \code{MxMatrix} 
objects and the \code{mxFitFunctionAlgebra} function. 
The \code{mxFitFunctionAlgebra} function must reference by name the \code{MxAlgebra} object to be evaluated.

The following operators and functions are supported in mxAlgebra:

Operators

\describe{
\item{\code{solve()}}{Inversion}
\item{\code{t()}}{Transposition}
\item{\code{^}}{Elementwise powering}
\item{\code{\%^\%}}{Kronecker powering}
\item{\code{+}}{Addition}
\item{\code{-}}{Subtraction}
\item{\code{\%*\%}}{Matrix Multiplication}
\item{\code{*}}{Elementwise product}
\item{\code{/}}{Elementwise division}
\item{\code{\%x\%}}{Kronecker product}
\item{\code{\%&\%}}{Quadratic product}
}

Functions

\describe{
\item{\code{cov2cor}}{Convert covariance matrix to correlation matrix}
\item{\code{chol}}{Cholesky Decomposition}
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
\item{\code{p2z}}{\emph{Standard}-normal quantile}
\item{\code{lgamma}}{Log-gamma function}
\item{\code{\link{eigenval}}}{Eigenvalues of a square matrix. Usage: eigenval(x); eigenvec(x); ieigenval(x); ieigenvec(x)}
\item{\code{\link{rvectorize}}}{Vectorize by row}
\item{\code{\link{cvectorize}}}{Vectorize by column}
\item{\code{\link{vech}}}{Half-vectorization}
\item{\code{\link{vechs}}}{Strict half-vectorization}
\item{\code{\link{vech2full}}}{Inverse half-vectorization}
\item{\code{\link{vechs2full}}}{Inverse strict half-vectorization}
\item{\code{\link{vec2diag}}}{Create matrix from a diagonal vector (similar to \link{diag}) }
\item{\code{\link{diag2vec}}}{Extract diagonal from matrix (similar to \link{diag}) }
\item{\code{\link{expm}}}{Matrix Exponential}
\item{\code{\link{logm}}}{Matrix Logarithm}
\item{\code{\link{omxExponential}}}{Matrix Exponential}
\item{\code{\link{omxMnor}}}{Multivariate Normal Integration}
\item{\code{\link{omxAllInt}}}{All cells Multivariate Normal Integration}
\item{\code{\link[=omxLogical]{omxNot}}}{Perform unary negation on a matrix}
\item{\code{\link[=omxLogical]{omxAnd}}}{Perform binary and on two matrices}
\item{\code{\link[=omxLogical]{omxOr}}}{Perform binary or on two matrices}
\item{\code{\link[=omxLogical]{omxGreaterThan}}}{Perform binary greater on two matrices}
\item{\code{\link[=omxLogical]{omxLessThan}}}{Perform binary less than on two matrices}
\item{\code{\link[=omxLogical]{omxApproxEquals}}}{Perform binary equals to (within a specified epsilon) on two matrices}
}}

\value{
    Returns a new \link{MxAlgebra} object.
}

\references{
The OpenMx User's guide can be found at http://openmx.psyc.virginia.edu/documentation.
}

\seealso{
\link{MxAlgebra} for the S4 class created by mxAlgebra. \link{mxFitFunctionAlgebra} for an objective function which takes an  MxAlgebra or MxMatrix object as the function to be minimized. 
\link{MxMatrix} and \link{mxMatrix} for objects which may be entered in the \code{expression} argument and the function that creates them. More information about the OpenMx package may be found \link[=OpenMx]{here}. 
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
model <- mxModel(model="AlgebraExample", A, B, C, D )
fit   <- mxRun(model)
mxEval(D, fit)


# Numbers in mxAlgebras are upgraded to 1x1 matrices
# Example of Kronecker powering (%^%) and multiplication (%*%)
A  <- mxMatrix(type="Full", nrow=3, ncol=3, value=c(1:9), name="A")
m1 <- mxModel(model="kron", A, mxAlgebra(A \%^\% 2, name="KroneckerPower"))
mxRun(m1)$KroneckerPower

# Running kron 
# mxAlgebra 'KroneckerPower' 
# $formula:  A %^% 2 
# $result:
#      [,1] [,2] [,3]
# [1,]    1   16   49
# [2,]    4   25   64
# [3,]    9   36   81

}
