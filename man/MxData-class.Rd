\name{MxData-class}
\alias{MxData}

\title{MxData Class}

\description{
   MxData is an S4 class. An MxData object is a \link[=Named-entity]{named entity}.
   New instances of this class can be created using the function \link{mxData}.
   MxData is an S4 class union. An MxData object is either \link{NULL} or a 
   MxNonNullData object.
}

\details{
The MxNonNullData class has the following slots:

  \tabular{rcl}{
	\tab \tab \cr
	name \tab - \tab The name of the object \cr
	observed \tab - \tab Either a matrix or a data frame \cr
	vector \tab - \tab A vector for means, or NA if missing \cr
	type \tab - \tab Either 'raw', 'cov', 'cor', or 'sscp' \cr
	numObs \tab - \tab The number of oberservations \cr
  }
  
The 'name' slot is the name of the MxData object. %Use of MxData objects in other functions in the \link{OpenMx} library may require reference by name.

The \sQuote{observed} slot is used to contain data, either as a matrix or as a data frame. Use of the data in this slot by other functions depends on the value of the 'type' slot. When 'type' is equal to 'cov', 'cor', or 'sscp', the data input into the 'matrix' slot should be a symmetric matrix or data frame.

The 'vector' slot is used to contain a vector of numeric values, which is used as a vector of means for MxData objects with 'type' equal to 'cov', 'cor', or 'sscp'. This slot may be used in estimation using the \link{mxMLObjective} function. 

The 'type' slot may take one of four supported values:
\itemize{
\item{raw}{The contents of the \sQuote{observed} slot are treated as raw data. Missing values are permitted and must be designated as the system missing value. The 'vector' and 'numObs' slots cannot be specified, as the 'vector' argument is not relevant and the 'numObs' argument is automatically populated with the number of rows in the data. Data of this type must use the \link{mxFIMLObjective} function as its objective function in MxModel objects, which deals with covariance estimation under full-information maximum likelihood.}

\item{cov}{The contents of the \sQuote{observed} slot are treated as a covariance matrix. The 'vector' argument is not required, but may be included for estimations involving means. The 'numObs' slot is required. Data of this type may use the \link{mxMLObjective}, or \link{mxRAMObjective} functions, depending on the specified model.}

\item{cor}{The contents of the \sQuote{observed} slot are treated as a correlation matrix. The 'vector' argument is not required, but may be included for estimations involving means. The 'numObs' slot is required.  Data of this type may use the \link{mxMLObjective}, or \link{mxRAMObjective} functions, depending on the specified model.}

\item{sscp}{The contents of the \sQuote{observed} slot are treated as a sums-of-squares and cross-products matrix. The 'vector' argument is not required, but may be included for estimations involving means. The 'numObs' slot is required.  Data of this type may use the \link{mxMLObjective}, or \link{mxRAMObjective} functions, depending on the specified model.}
}

The 'numObs' slot describes the number of observations in the data. If 'type' equals 'raw', then 'numObs' is automatically populated as the number of rows in the matrix or data frame in the \sQuote{observed} slot. If 'type' equals 'cov', 'cor', or 'sscp', then this slot must be input using the 'numObs' argument in the \link{mxData} function when the MxData argument is created.

MxData objects may not be included in \link{MxAlgebra} objects or use the \link{mxAlgebraObjective} function. If these capabilities are desired, data should be appropriately input or transformed using the \link{mxMatrix} and \link{mxAlgebra} functions.

While column names are stored in the \sQuote{observed} slot of MxData objects, these names are not recognized as variable names in \link{MxPath} objects. Variable names must be specified using the 'manifestVars' argument of the \link{mxModel} function prior to use in MxPath objects.

The mxData function does not currently place restrictions on the size, shape, or symmetry of matrices input into the \sQuote{observed} argument. While it is possible to specify MxData objects as covariance, correlation or sscp matrices that do not have the properties commonly associated with these matrices, failure to correctly specify these matrices will likely lead to problems in model estimation.
}

\references{
The Mx manual may be found at http://www.vcu.edu/mx/documentation.html.
}

\seealso{
\link{mxData} for creating MxData objects, \link{matrix} and \link{data.frame} for objects which may be entered as arguments in the 'matrix' slot. More information about the OpenMx package may be found \link[=OpenMx]{here}.
}
