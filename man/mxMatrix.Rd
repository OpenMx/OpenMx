\name{mxMatrix}
\alias{mxMatrix}

\title{Function To Create MxMatrix Object}

\description{
  This functions creates a new \link{MxMatrix} object.
}

\usage{
mxMatrix(type = "Full", values = NA, free = FALSE, labels = NA,
	lbound = NA, ubound = NA, nrow = NA, ncol = NA, 
	byrow = FALSE, name = NA)
}

\arguments{
  \item{type}{a character string indicating the matrix type, where type indicates the range of values and equalities in the matrix. Must be one of:  'Diag', 'Full', 'Iden', 'Symm', 'Unit', or 'Zero'.}
  \item{values}{a vector or matrix of numeric starting values. By default, all values are set to zero.}
  \item{free}{a vector or matrix of logicals for free parameter specification. A single 'TRUE' or 'FALSE' will set all allowable variables to free or fixed, respectively.}
  \item{labels}{a vector or matrix of characters for variable label specification.} 
  \item{lbound}{a vector or matrix of numeric lower bounds. Default bounds are specified with an NA.}
  \item{ubound}{a vector or matrix of numeric upper bounds. Default bounds are specified with an NA.}  
  \item{nrow}{the desired number of rows. One of 'nrow' and 'ncol' is required when 'values', 'free', and 'labels' arguments are not matrices.}
  \item{ncol}{the desired number of columns. One of 'nrow' and 'ncol' is required when 'values', 'free', and 'labels' arguments are not matrices.}
  \item{byrow}{logical. If 'FALSE' (default), the 'values', 'free', and 'labels' matrices are populated by column rather than by row. }
  \item{name}{an optional character string indicating the name of the MxMatrix object created by the mxModel function.}
}

\details{
The mxMatrix function creates \link{MxMatrix} objects, which consist of a pair of matrices and a 'type' argument. The 'values' matrix is made up of numeric elements whose usage and capabilities in other functions are defined by the 'free' matrix. If an element is specified as a fixed parameter in the 'free' matrix, then the element in the 'values' matrix is treated as a constant value and cannot be altered or updated by an objective function when included in an \link{mxRun} function. If an element is specified as a free parameter in the 'free' matrix, the element in the 'value' matrix is considered a starting value and can be changed by an objective function when included in an \link{mxRun} function. Free parameters are specified with a character string, non-zero numeric value, or 'NA'; fixed parameters are specified with a numeric zero.

Objects created by the mxMatrix function are of a specific 'type', which specifies the number and location of parameters in the 'spec' matrix and the possible values in the 'value' matrix. Input 'value', 'free', and 'label' matrices must be of appropriate shape and have appropriate values for the matrix type requested. Six types of matrices are supported:

\itemize{
\item{'Diag'}{matrices must be square, and only elements on the principle diagonal may be specified as free parameters or take non-zero values. All other elements are required to be fixed parameters with a value of 0.}

\item{'Full'}{matrices may be either rectangular or square, and all elements in the matrix may be freely estimated. This type is the default for the mxMatrix() function.}

\item{'Iden'}{matrices must be square, and consist of no free parameters. Matrices of this type have a value of 1 for all entries on the principle diagonal and the value 0 in all off-diagonal entries.}

\item{'Symm'}{matrices must be square, and elements in the principle diagonal and lower triangular portion of the matrix may be free parameters of any value. Elements in the upper triangular portion of the matrix are constrained to be equal to those in the lower triangular portion, such that the value and parameter specificiation of the element in row i and column j is identical to to the value and specification of the element in row j and column i. Values and parameters specified in one element of the matrix are automatically applied to equivalent elements.}

\item{'Unit'}{matrices may be either rectangular or square, and contain no free parameters. All elements in matrices of this type have a value of 1 for all elements.}

\item{'Zero'}{matrices may be either rectangular or square, and contain no free parameters. All elements in matrices of this type have a value of 0 for all elements.}
}

}

\value{
  Returns a new \link{MxMatrix} object, which consists of a 'value' matrix of numeric starting values, a 'free' matrix describing free parameter specification, and a 'labels' matrix of labels for the variable names. This \link{MxMatrix} object can be used as an argument in the \link{mxAlgebra}, \link{mxBounds}, \link{mxConstraint} and \link{mxModel} functions.
}

\references{
More information about these matrix types can be found in Chapter 4 of the Mx manual, which can be found at http://www.vcu.edu/mx/documentation.html.
}

\seealso{
\link{MxMatrix} for the S4 class created by mxMatrix. More information about the OpenMx package may be found \link[=OpenMx]{here}. 
}

\examples{
#Create a 3x3 identity matrix
i <- mxMatrix(type="Iden", nrow=3, ncol=3)

#Create a full 4x2 matrix from existing value and specification matrices with all free parameters
value <- matrix(c(1,2,3,4,5,6,7,8), nrow=4)
f <- mxMatrix(value=value, free = TRUE)

#Create a 3x3 symmetric matrix with free off-diagonal parameters and starting values
s <- mxMatrix(type="Symm", 
	c(1, .8, .8, .8, 1, .8, .8, .8, 1),
	c(FALSE, TRUE, TRUE, TRUE, FALSE, 
	TRUE, TRUE, TRUE, FALSE), 
	c(NA, "free1", "free2", "free1", NA, 
	"free3", "free2", "free3", NA), nrow=3)
}
