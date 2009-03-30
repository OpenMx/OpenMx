\name{mxMatrix}
\alias{mxMatrix}

\title{Function To Create MxMatrix Object}

\description{
  This functions creates a new \link{MxMatrix} object.
}

\usage{
mxMatrix(type = "Full", values = NA, specification = NA,
         nrow = NA, ncol = NA, byrow = FALSE,
         free = FALSE, name = NA}
}

\arguments{
  \item{type}{a character string indicating the matrix type, where type indicates the range of values and equalities in the matrix. Must be one of:  'Diag', 'Full', 'Iden', 'Symm', 'Unit', or 'Zero'.}
  \item{values}{a vector or matrix of numeric starting values. By default, all values are set to zero.}
  \item{specification}{a vector or matrix of character strings used for free parameter specification. Defaults can be set by the 'free' argument as either all fixed or all free.}
  \item{nrow}{the desired number of rows. Required when 'values' and 'specification' arguments are not matrices. Must be equal to the 'ncol' argument when matrix is of type 'Diag', 'Iden' and 'Symm'.}
  \item{ncol}{the desired number of columns. Required when 'values' and 'specification' arguments are not matrices. Must be equal to the 'nrow' argument when matrix is of type 'Diag', 'Iden' and 'Symm'.}
  \item{byrow}{logical. If 'FALSE' (default), the 'values' and 'specification' matrices are populated by row rather than by column, when 'values' and 'specification' arguments are not matrices. }
  \item{free}{logical. If 'FALSE' (default), the 'specification' matrix is populated with fixed parameters, represented by a zero. If 'TRUE', the 'specification' matrix is populated with free parameters. Must be 'FALSE' if 'type' is equal to 'Iden', 'Unit', or 'Zero'.}
  \item{name}{an optional character string indicating the name of the MxMatrix object created by the mxModel function.}
}

\details{
\link{MxMatrix} objects consist of a pair of matrices and a 'type' argument. The 'values' matrix is made up of numeric elements whose usage and capabilities in other functions are defined by the 'specification' matrix. If an element is specified as a fixed parameter in the 'specification' matrix, then the element in the 'values' matrix is treated as a constant value and cannot be altered or updated by an objective function when included in an (mxJob) function. If an element is specified as a free parameter in the 'specification' matrix, the element in the 'value' matrix is considered a starting value and can be changed by an objective function when included in an (mxJob) function.

Objects created by the mxMatrix function are of a specific 'type', which specifies the number and location of free parameters in the 'specification' matrix and the possible values in the 'value' matrix. Six types of matrices are supported:

\item{'Diag'}{matrices must be square, and only elements on the principle diagonal may be specified as free parameters or take non-zero values. All other elements are required to be fixed parameters with a value of 0.}

\item{'Full'}{matrices may be either rectangular or square, and all elements in the matrix may be freely estimated. This type is the default for the mxMatrix() function.}

\item{'Iden'}{matrices must be square, and consist of no free parameters. Matrices of this type have a value of 1 for all entries on the principle diagonal and the value 0 in all off-diagonal entries.}

\item{'Symm'}{matrices must be square, and elements in the principle diagonal and lower triangular portion of the matrix may be free parameters of any value. Elements in the upper triangular portion of the matrix are constrained to be equal to those in the lower triangular portion, such that the value and parameter specificiation of the element in row i and column j is identical to to the value and specification of the element in row j and column i. Values and parameters specified in one element of the matrix are automatically applied to equivalent elements.}

\item{'Unit'}{matrices may be either rectangular or square, and contain no free parameters. All elements in matrices of this type have a value of 1 for all elements.}

\item{'Zero'}{matrices may be either rectangular or square, and contain no free parameters. All elements in matrices of this type have a value of 0 for all elements.}

}

\value{
  Returns a new \link{MxMatrix} object, which consists of a 'value' matrix of numeric starting values, a 'specification' matrix describing free parameter specification, and a 'type' argument, which must be either "Diag", "Full", "Iden", "Symm", "Unit",  or "Zero". This \link{MxMatrix} object can be used as an argument in the \link{mxAlgebra}, \link{mxBounds}, \link{mxConstraint} and \link{mxModel} functions.\\
}

\references{
More information about these matrix types can be found in Chapter 4 of the Mx manual, which can be found at http://www.vcu.edu/mx/documentation.html.
}

\examples{
#Create a 3x3 identity matrix
i<-mxMatrix(type="Iden",nrow=3,ncol=3)
}