\name{mxMatrix}
\alias{mxMatrix}

\title{Function To Create MxMatrix Object}

\description{
   This function creates a new \link{MxMatrix} object.
}

\usage{
mxMatrix(type = "Full", values = NA, specification = NA,
         nrow = NA, ncol = NA, byrow = FALSE, 
         free = FALSE, name = NA)
}

\arguments{
   \item{type}{A character string indicating the matrix type. 
		Must be either "Diag", "Full", "Iden", "Symm", "Unit", or "Zero".}
   \item{values}{A vector or matrix of starting values.}
   \item{specification}{A vector or matrix for free parameter specifications.}
   \item{nrow}{The desired number of rows, 
		when 'values' and 'specification' are not matrices.}
   \item{ncol}{The desired number of cols, 
		when 'values' and 'specification' are not matrices.}
   \item{byrow}{logical.  If 'FALSE' (the default) the matrices are filled by
		columns, otherwise the matrices are filled by rows.}
   \item{free}{logical.  If 'FALSE' (the default) the matrix is fixed, otherwise
		the matrix is free.}
   \item{name}{A character string indicating the name of the object. Optional.}
}

\value{
	Returns a new \link{MxMatrix} object.
}

\examples{

   # create a 2x3 full matrix with free parameters
   A <- mxMatrix("Full", values = 0, nrow = 2, ncol = 3, 
		free = TRUE, name = "A")

   # create a 3x3 fixed identity matrix
   I <- mxMatrix("Iden", nrow = 3, name = "I")

   # create a 3x3 symmetric matrix with no free parameters, fill by row
   C <- mxMatrix("Symm", values = 1:6, nrow = 3, 
		byrow = TRUE, name = "C")
}
