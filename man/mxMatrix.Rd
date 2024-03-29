%
%   Copyright 2007-2021 by the individuals mentioned in the source code history
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

\name{mxMatrix}
\alias{mxMatrix}

\title{Create MxMatrix Object}

\description{
  This function creates a new \link{MxMatrix} object.
}

\usage{
mxMatrix(type = c('Full', 'Diag', 'Iden', 'Lower',
	'Sdiag', 'Stand', 'Symm', 'Unit', 'Zero'), nrow = NA, ncol = NA, 
    free = FALSE, values = NA, labels = NA, lbound = NA, 
    ubound = NA, byrow = getOption('mxByrow'), dimnames = NA, name = NA, 
    condenseSlots=getOption('mxCondenseMatrixSlots'),
    ..., joinKey=as.character(NA), joinModel=as.character(NA))
}

\arguments{
  \item{type}{A character string indicating the matrix type, where type indicates the range of values and equalities in the matrix. Must be one of:  \sQuote{Diag}, \sQuote{Full}, \sQuote{Iden}, \sQuote{Lower}, \sQuote{Sdiag}, \sQuote{Stand}, \sQuote{Symm}, \sQuote{Unit}, or \sQuote{Zero}.}
  \item{nrow}{Integer; the desired number of rows. One or both of \sQuote{nrow} and \sQuote{ncol} is required when \sQuote{values}, \sQuote{free}, \sQuote{labels}, \sQuote{lbound}, and \sQuote{ubound} arguments are not matrices, depending on the desired \link{MxMatrix} type.}
  \item{ncol}{Integer; the desired number of columns. One or both of \sQuote{nrow} and \sQuote{ncol} is required when \sQuote{values}, \sQuote{free}, \sQuote{labels}, \sQuote{lbound}, and \sQuote{ubound} arguments are not matrices, depending on the desired \link{MxMatrix} type.}
  \item{free}{A vector or matrix of logicals for free parameter specification. A single \sQuote{TRUE} or \sQuote{FALSE} will set all allowable variables to free or fixed, respectively.}
  \item{values}{A vector or matrix of numeric starting values. By default, all values are set to zero.}
  \item{labels}{A vector or matrix of characters for variable label specification.} 
  \item{lbound}{A vector or matrix of numeric lower bounds. Default bounds are specified with an \code{NA}.}
  \item{ubound}{A vector or matrix of numeric upper bounds. Default bounds are specified with an \code{NA}.}  
  \item{byrow}{Logical; defaults to value of global \link[base:options]{option} 'mxByRow'. If \code{FALSE} (default), the \sQuote{values}, \sQuote{free}, \sQuote{labels}, \sQuote{lbound}, and \sQuote{ubound} matrices are populated by column rather than by row. }
  \item{dimnames}{List. The dimnames attribute for the matrix: a list of length 2 giving the row and column names respectively. An empty list is treated as NULL, and a list of length one as row names. The list can be named, and the list names will be used as names for the dimensions.}
  \item{name}{An optional character string indicating the name of the \link{MxMatrix} object.}
  \item{condenseSlots}{Logical; defaults to value of global \link[base:options]{option} 'mxCondenseMatrixSlots'.  If \code{TRUE}, then the resulting \link{MxMatrix} will "condense" its \sQuote{labels}, \sQuote{free}, \sQuote{lbound}, and \sQuote{ubound} down to 1x1 matrices if they contain only \code{FALSE} (\sQuote{free}) or \code{NA} (the other three).  If \code{FALSE}, those four matrices and the \sQuote{values} matrix will all be of equal dimensions.   }
\item{...}{Not used.  Forces remaining arguments to be specified by
  name.}
\item{joinKey}{The name of the column in current model's raw data that
  is used as a foreign key to match against the primary key in the
  joinModel's raw data.}
\item{joinModel}{The name of the model that this matrix joins against.}
}

\details{
The mxMatrix function creates \link{MxMatrix} objects, which consist of five matrices and a \sQuote{type} argument. The \sQuote{values} matrix is made up of numeric elements whose usage and capabilities in other functions are defined by the \sQuote{free} matrix. If an element is specified as a fixed parameter in the \sQuote{free} matrix, then the element in the \sQuote{values} matrix is treated as a constant value and cannot be altered or updated by an objective function when included in an \link{mxRun} function. If an element is specified as a free parameter in the \sQuote{free} matrix, the element in the \sQuote{value} matrix is considered a starting value and can be changed by an objective function when included in an \link{mxRun} function. %Free parameters are specified with a character string, non-zero numeric value, or \sQuote{NA}; fixed parameters are specified with a numeric zero.

Element labels beginning with \code{'data.'} can be used if the \link{MxMatrix} is to be used in an \link{MxModel} object that has a raw dataset (i.e., an \link{MxData} object of \code{type="raw"}).  Such a label instructs OpenMx to use a particular column of the raw dataset to fill in the value of that element.  For historical reasons, the variable contained in that column is called a "definition variable."  For example, if an \link{MxMatrix} element has the label \code{'data.x'}, then OpenMx will use the first value of the data column named "x" when evaluating the fitfunction for the first row, and will use the second value of column "x" when evaluating the fitfunction for the second row, and so on.  After the call to \code{\link{mxRun}()}, the values for elements labeled with \code{'data.x'} are returned as the value from the first (i.e., first before any automated sorting is done) element of column "x" in the data. 

Objects created by the \code{mxMatrix()} function are of a specific \sQuote{type}, which specifies the number and location of parameters in the \sQuote{labels} matrix and the starting values in the \sQuote{values} matrix. Input \sQuote{values}, \sQuote{free}, and \sQuote{labels} matrices must be of appropriate shape and have appropriate values for the matrix type requested. Nine types of matrices are supported:

\tabular{ll}{
\sQuote{Diag}   \tab matrices must be square, and only elements on the principal diagonal may be specified as free parameters or take non-zero values. All other elements are required to be fixed parameters with a value of 0. \cr

\sQuote{Full}   \tab matrices may be either rectangular or square, and all elements in the matrix may be freely estimated. This type is the default for the mxMatrix() function. \cr

\sQuote{Iden}   \tab matrices must be square, and consist of no free parameters. Matrices of this type have a value of 1 for all entries on the principal diagonal and the value 0 in all off-diagonal entries. \cr

\sQuote{Lower}  \tab matrices must be square, with a value of 0 for all entries in the upper triangle and no free parameters in the upper triangle. \cr

\sQuote{Sdiag}  \tab matrices must be square, with a value of 0 for all entries in the upper triangle and along the diagonal. No free parameters in the upper triangle or along the diagonal. \cr

\sQuote{Symm}   \tab matrices must be square, and elements in the principle diagonal and lower triangular portion of the matrix may be free parameters of any value. Elements in the upper triangular portion of the matrix are constrained to be equal to those in the lower triangular portion, such that the value and parameter specification of the element in row i and column j is identical to to the value and specification of the element in row j and column i. \cr

\sQuote{Stand}  \tab matrices are symmetric matrices (see 'Symm') with 1's along the main diagonal. \cr

\sQuote{Unit}   \tab matrices may be either rectangular or square, and contain no free parameters. All elements in matrices of this type have a value of 1 for all elements. \cr

\sQuote{Zero}   \tab matrices may be either rectangular or square, and contain no free parameters. All elements in matrices of this type have a value of 0 for all elements. \cr
}

When \sQuote{type} is \sQuote{Lower} or \sQuote{Symm}, then the arguments to \sQuote{free}, \sQuote{values}, \sQuote{labels}, \sQuote{lbound}, or \sQuote{ubound} may be vectors of length \eqn{N * (N + 1) / 2}, where N is the number of rows and columns of the matrix. When \sQuote{type} is \sQuote{Sdiag} or \sQuote{Stand}, then the arguments to \sQuote{free}, \sQuote{values}, \sQuote{labels}, \sQuote{lbound}, or \sQuote{ubound} may be vectors of length \eqn{N * (N - 1) / 2}. 

}

\value{
  Returns a new \link{MxMatrix} object, which consists of a \sQuote{values} matrix of numeric starting values, a \sQuote{free} matrix describing free parameter specification, a \sQuote{labels} matrix of labels for the variable names, and \sQuote{lbound} and \sQuote{ubound} matrices of the lower and upper parameter bounds. This \link{MxMatrix} object can be used as an argument in the \code{\link{mxAlgebra}()}, \code{\link{mxBounds}()}, \code{\link{mxConstraint}()} and \code{\link{mxModel}()} functions.
}

\references{
The OpenMx User's guide can be found at \url{https://openmx.ssri.psu.edu/documentation/}.
}

\seealso{
\link{MxMatrix} for the S4 class created by mxMatrix. More information about the OpenMx package may be found \link[=OpenMx]{here}. 
}

\examples{
# Create a 3 x 3 identity matrix

idenMatrix <- mxMatrix(type = "Iden", nrow = 3, 
    ncol = 3, name = "I")

# Create a full 4 x 2 matrix from existing 
# value matrix with all free parameters

vals <- matrix(1:8, nrow = 4)
fullMatrix <- mxMatrix(type = "Full", values = vals, 
    free = TRUE, name = "foo")
 
# Create a 3 x 3 symmetric matrix with free off-
# diagonal parameters and starting values

symmMatrix <- mxMatrix(type = "Symm", nrow = 3, ncol = 3,
    free = c(FALSE, TRUE, TRUE, FALSE, TRUE, FALSE),
    values = c(1, .8, .8, 1, .8, 1),
    labels = c(NA, "free1", "free2", NA, "free3", NA),
    name = "bar")

# Create an mxMatrix from a character matrix.  All numbers are
# interpreted as fixed and non-numbers are interpreted as free
# parameters.

matrixFromChar <- function(inputm, name=NA) {
  inputmFixed <- suppressWarnings(matrix(
    as.numeric(inputm),nrow = nrow(inputm), ncol = ncol(inputm)))
  inputmCharacter <- inputm
  inputmCharacter[!is.na(inputmFixed)] <- NA
  mxMatrix(nrow=nrow(inputm), ncol=ncol(inputm),
           free=!is.na(inputmCharacter),
           values=inputmFixed,
           labels=inputmCharacter,
           dimnames=dimnames(inputm), name=name)
}


# Demonstrate some of the behavior of the condensed slots
# Create a 3x3 matrix with condensed slots

a <- mxMatrix('Full', 3, 3, values=1, condenseSlots=TRUE)
a@free # at operator returns the stored 1x1 matrix
a$free # dollar operator constructs full matrix for printing

# assignment with the dollar operator
#  de-condenses the slots to create the
#  full 3x3 matrix
a$free[1,1] <- TRUE
a@free



}
