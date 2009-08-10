\name{mxMatrix}
\alias{mxMatrix}

\title{Function To Create MxMatrix Object}

\description{
  This functions creates a new \link{MxMatrix} object.
}

\usage{
mxMatrix(type = "Full", nrow = NA, ncol = NA, 
	free = FALSE, values = NA, labels = NA, lbound = NA, 
	ubound = NA, byrow = FALSE, dimnames = NA, name = NA)
}

\arguments{
  \item{type}{a character string indicating the matrix type, where type indicates the range of values and equalities in the matrix. Must be one of:  \sQuote{Diag}, \sQuote{Full}, \sQuote{Iden}, \sQuote{Lower}, \sQuote{Sdiag}, \sQuote{Stand}, \sQuote{Symm}, \sQuote{Unit}, or \sQuote{Zero}.}
  \item{nrow}{the desired number of rows. One or both of \sQuote{nrow} and \sQuote{ncol} is required when \sQuote{values}, \sQuote{free}, \sQuote{labels}, \sQuote{lbound}, and \sQuote{ubound} arguments are not matrices, depending on the matrix type.}
  \item{ncol}{the desired number of columns. One or both of \sQuote{nrow} and \sQuote{ncol} is required when \sQuote{values}, \sQuote{free}, \sQuote{labels}, \sQuote{lbound}, and \sQuote{ubound} arguments are not matrices, depending on the matrix type.}
  \item{free}{a vector or matrix of logicals for free parameter specification. A single \sQuote{TRUE} or \sQuote{FALSE} will set all allowable variables to free or fixed, respectively.}
  \item{values}{a vector or matrix of numeric starting values. By default, all values are set to zero.}
  \item{labels}{a vector or matrix of characters for variable label specification.} 
  \item{lbound}{a vector or matrix of numeric lower bounds. Default bounds are specified with an NA.}
  \item{ubound}{a vector or matrix of numeric upper bounds. Default bounds are specified with an NA.}  
  \item{byrow}{logical. If \sQuote{FALSE} (default), the \sQuote{values}, \sQuote{free}, \sQuote{labels}, \sQuote{lbound}, and \sQuote{ubound} matrices are populated by column rather than by row. }
  \item{dimnames}{list. The dimnames attribute for the matrix: a list of length 2 giving the row and column names respectively. An empty list is treated as NULL, and a list of length one as row names. The list can be named, and the list names will be used as names for the dimensions.}
  \item{name}{an optional character string indicating the name of the MxMatrix object created by the mxModel function.}
}

\details{
The mxMatrix function creates \link{MxMatrix} objects, which consist of a pair of matrices and a \sQuote{type} argument. The \sQuote{values} matrix is made up of numeric elements whose usage and capabilities in other functions are defined by the \sQuote{free} matrix. If an element is specified as a fixed parameter in the \sQuote{free} matrix, then the element in the \sQuote{values} matrix is treated as a constant value and cannot be altered or updated by an objective function when included in an \link{mxRun} function. If an element is specified as a free parameter in the \sQuote{free} matrix, the element in the \sQuote{value} matrix is considered a starting value and can be changed by an objective function when included in an \link{mxRun} function. Free parameters are specified with a character string, non-zero numeric value, or \sQuote{NA}; fixed parameters are specified with a numeric zero.

Objects created by the mxMatrix function are of a specific \sQuote{type}, which specifies the number and location of parameters in the \sQuote{labels} matrix and the starting values in the \sQuote{values} matrix. Input \sQuote{values}, \sQuote{free}, and \sQuote{labels} matrices must be of appropriate shape and have appropriate values for the matrix type requested. Nine types of matrices are supported:

\tabular{ll}{
\sQuote{Diag}   \tab matrices must be square, and only elements on the principle diagonal may be specified as free parameters or take non-zero values. All other elements are required to be fixed parameters with a value of 0. \cr

\sQuote{Full}   \tab matrices may be either rectangular or square, and all elements in the matrix may be freely estimated. This type is the default for the mxMatrix() function. \cr

\sQuote{Iden}   \tab matrices must be square, and consist of no free parameters. Matrices of this type have a value of 1 for all entries on the principle diagonal and the value 0 in all off-diagonal entries. \cr

\sQuote{Lower}  \tab matrices must be square, with a value of 0 for all entries in the upper triangle and no free parameters in the upper triangle. \cr

\sQuote{Sdiag}  \tab matrices must be square, with a value of 0 for all entries in the upper triangle and along the diagonal. No free parameters in the upper triangle or along the diagonal. \cr

\sQuote{Symm}   \tab matrices must be square, and elements in the principle diagonal and lower triangular portion of the matrix may be free parameters of any value. Elements in the upper triangular portion of the matrix are constrained to be equal to those in the lower triangular portion, such that the value and parameter specificiation of the element in row i and column j is identical to to the value and specification of the element in row j and column i. Values and parameters specified in one element of the matrix are automatically applied to equivalent elements. \cr

\sQuote{Stand}  \tab matrices are symmetric matrices (see 'Symm') with 1's along the main diagonal. \cr

\sQuote{Unit}   \tab matrices may be either rectangular or square, and contain no free parameters. All elements in matrices of this type have a value of 1 for all elements. \cr

\sQuote{Zero}   \tab matrices may be either rectangular or square, and contain no free parameters. All elements in matrices of this type have a value of 0 for all elements. \cr
}

}

\value{
  Returns a new \link{MxMatrix} object, which consists of a \sQuote{values} matrix of numeric starting values, a \sQuote{free} matrix describing free parameter specification, a \sQuote{labels} matrix of labels for the variable names, and \sQuote{lbound} and \sQuote{ubound} matrices of the lower and upper parameter bounds. This \link{MxMatrix} object can be used as an argument in the \link{mxAlgebra}, \link{mxBounds}, \link{mxConstraint} and \link{mxModel} functions.
}

\references{
The OpenMx User's guide can be found at http://openmx.psyc.virginia.edu/documentation.
}

\seealso{
\link{MxMatrix} for the S4 class created by mxMatrix. More information about the OpenMx package may be found \link[=OpenMx]{here}. 
}

\examples{
#Create a 3x3 identity matrix
i <- mxMatrix(type="Iden", 3, 3)

#Create a full 4x2 matrix from existing value and specification matrices with all free parameters
value <- matrix(c(1,2,3,4,5,6,7,8), nrow=4)
f <- mxMatrix(value=value, free = TRUE)

#Create a 3x3 symmetric matrix with free off-diagonal parameters and starting values
s <- mxMatrix(type="Symm", 3, 3,
	c(FALSE, TRUE, TRUE, TRUE, FALSE, 
	TRUE, TRUE, TRUE, FALSE),
	c(1, .8, .8, .8, 1, .8, .8, .8, 1),
	c(NA, "free1", "free2", "free1", NA, 
	"free3", "free2", "free3", NA))
}
