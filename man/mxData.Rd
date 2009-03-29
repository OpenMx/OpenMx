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
   \item{matrix}{This argument is either a matrix or a data.frame}
   \item{type}{This argument is a character string.  It must be one of 'raw', 'cov', 'cor', or 'sscp'.}
   \item{vector}{A means vector. Optional.}
   \item{numObs}{The number of observations.  If type is not 'raw', then this field is required.}
   \item{name}{A character string indicating the name of the object. Optional.}
}

\value{
	Returns a new \link{MxData} object.
}

\examples{
	
# Create an anonymous model
model <- mxModel()

covMatrix <- matrix( c(0.77642931, 0.39590663, 
    0.39590663, 0.49115615), 
    nrow = 2, ncol = 2, byrow = TRUE)

data <- mxData(covMatrix, 'cov', numObs = 100, 
    name = 'covariance matrix')

model <- mxModel(model, data) 

}
