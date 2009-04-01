\name{mxData}
\alias{mxData}

\title{Function To Create MxData Object}

\description{
   This function creates a new \link{MxData} object.
}

\usage{
mxData(matrix, type, vector = NA, numObs = NA, name = NA)
}

\arguments{
   \item{matrix}{A matrix or data.frame which provides data to the MxData object.}
   \item{type}{A character string defining the type of data in the 'matrix' argument. Must be one of 'raw', 'cov', 'cor', or 'sscp'.}
   \item{vector}{An optional vector of means for use when 'type' is 'cov', 'cor' or 'sscp'.}
   \item{numObs}{The number of observations in the data supplied in the 'matrix' argument. Required unless 'type' equals 'raw'.}
   \item{name}{An optional character string indicating the name of the resulting MxData object.}
}

\details{
The mxData function creates \link{MxData} objects, which can be used as arguments in \link{MxModel} objects. The 'matrix' argument may take either a data.frame or a matrix, which is then described with the 'type' argument. Data types describe compatibility and usage with objective functions in MxModel objects. Four different data types are supported:

\item{raw}{The contents of the 'matrix' argument are treated as raw data. Missing values are permitted and must be designated as the system missing value. The 'vector' and 'numObs' arguments cannot be specified, as the 'vector' argument is not relevant and the 'numObs' argument is automatically populated with the number of rows in the data. Data of this type must use the \link{mxFIMLObjective} function as its objective function in MxModel objects, which deals with covariance estimation under full-information maximum likelihood.}

\item{cov}{The contents of the 'matrix' argument are treated as a covariance matrix. The 'vector' argument is not required, but may be included for estimations involving means. The 'numObs' argument is required, which should reflect the number of observations or rows in the data described by the covariance matrix. Data of this type may use the \link{mxMLObjective}, or \link{mxRAMObjective} functions, depending on the specified model.}

\item{cor}{The contents of the 'matrix' argument are treated as a correlation matrix. The 'vector' argument is not required, but may be included for estimations involving means. The 'numObs' argument is required, which should reflect the number of observations or rows in the data described by the covariance matrix. Data of this type may use the \link{mxMLObjective}, or \link{mxRAMObjective} functions, depending on the specified model.}

\item{sscp}{The contents of the 'matrix' argument are treated as a sums-of-squares and cross-products matrix. The 'vector' argument is not required, but may be included for estimations involving means. The 'numObs' argument is required, which should reflect the number of observations or rows in the data described by the covariance matrix. Data of this type may use the \link{mxMLObjective}, or \link{mxRAMObjective} functions, depending on the specified model.}

MxData objects may not be included in \link{MxAlgebra} objects or use the \link{mxAlgebraObjective} function. If these capabilities are desired, data should be appropriately input or transformed using the \link{mxMatrix} and \link{mxAlgebra} functions.

While column names are stored in the 'matrix' slot of MxData objects, these names are not recognized as variable names in \link{MxPath} objects. Variable names must be specified using the 'manifestVars' argument of the \link{mxModel} function prior to use in MxPath objects.

The mxData function does not currently place restrictions on the size, shape, or symmetry of matrices input into the 'matrix' argument. While it is possible to specify MxData objects as covariance, correlation or sscp matrices that do not have the properties commonly associated with these matrices, failure to correctly specify these matrices will likely lead to problems in model estimation.
}

\value{
	Returns a new \link{MxData} object.
}

\examples{	
#Create a covariance matrix
covMatrix <- matrix( c(0.77642931, 0.39590663, 
    0.39590663, 0.49115615), 
    nrow = 2, ncol = 2, byrow = TRUE)

#Create an MxData object including that covariance matrix
data <- mxData(covMatrix, 'cov', numObs = 100, 
    name = 'covmatrix')

model <- mxModel(model, data) 

}
