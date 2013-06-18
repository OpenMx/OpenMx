%
%   Copyright 2007-2013 The OpenMx Project
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

\name{mxData}
\alias{mxData}

\title{Create MxData Object}

\description{
   This function creates a new \link{MxData} object.
}

\usage{
mxData(observed, type, means = NA, numObs = NA)
}

\arguments{
   \item{observed}{A matrix or data.frame which provides data to the MxData object.}
   \item{type}{A character string defining the type of data in the \sQuote{observed} argument. Must be one of \dQuote{raw}, \dQuote{cov}, or \dQuote{cor}.}
   \item{means}{An optional vector of means for use when \sQuote{type} is \dQuote{cov}, or \dQuote{cor}.}
   \item{numObs}{The number of observations in the data supplied in the \sQuote{observed} argument. Required unless \sQuote{type} equals \dQuote{raw}.}
}

\details{
The mxData function creates \link{MxData} objects, which can be used as arguments in \link{MxModel} objects. The \sQuote{observed} argument may take either a data frame or a matrix, which is then described with the \sQuote{type} argument. Data types describe compatibility and usage with objective functions in MxModel objects. Four different data types are supported:
\describe{
\item{raw}{The contents of the \sQuote{observed} argument are treated as raw data. Missing values are permitted and must be designated as the system missing value. The \sQuote{means} and \sQuote{numObs} arguments cannot be specified, as the \sQuote{means} argument is not relevant and the \sQuote{numObs} argument is automatically populated with the number of rows in the data. Data of this type must use the \link{mxFIMLObjective} function as its objective function in MxModel objects, which deals with covariance estimation under full-information maximum likelihood.}

\item{cov}{The contents of the \sQuote{observed} argument are treated as a covariance matrix. The \sQuote{means} argument is not required, but may be included for estimations involving means. The \sQuote{numObs} argument is required, which should reflect the number of observations or rows in the data described by the covariance matrix. Data of this type may use the \link{mxMLObjective}, or \link{mxRAMObjective} functions, depending on the specified model.}

\item{cor}{The contents of the \sQuote{observed} argument are treated as a correlation matrix. The \sQuote{means} argument is not required, but may be included for estimations involving means. The \sQuote{numObs} argument is required, which should reflect the number of observations or rows in the data described by the covariance matrix. Data of this type may use the \link{mxMLObjective}, or \link{mxRAMObjective} functions, depending on the specified model.}
}

MxData objects may not be included in \link{MxAlgebra} objects or use the \link{mxAlgebraObjective} function. If these capabilities are desired, data should be appropriately input or transformed using the \link{mxMatrix} and \link{mxAlgebra} functions.

While column names are stored in the \sQuote{observed} slot of MxData objects, these names are not recognized as variable names in \link{MxPath} objects. Variable names must be specified using the \sQuote{manifestVars} argument of the \link{mxModel} function prior to use in MxPath objects.

The mxData function does not currently place restrictions on the size, shape, or symmetry of matrices input into the \sQuote{observed} argument. While it is possible to specify MxData objects as covariance or correlation  matrices that do not have the properties commonly associated with these matrices, failure to correctly specify these matrices will likely lead to problems in model estimation.

OpenMx uses the names of variables to map them onto the objective functions and other elements associated with your model. For data.frames, ensure you have set the names(). For matrices set names using, for instance, row.names=c("your", "columns"). Covariance  and cor matrices need to have both the row and column names set and these must be identical, for instance by using dimnames=list(varNames, varNames).


}

\value{
	Returns a new \link{MxData} object.
}

\references{
The OpenMx User's guide can be found at http://openmx.psyc.virginia.edu/documentation.
}

\seealso{
\link{MxData} for the S4 class created by mxData. \link{matrix} and \link{data.frame} for objects which may be entered as arguments in the \sQuote{observed} slot. More information about the OpenMx package may be found \link[=OpenMx]{here}. 
}

\examples{	
#Create a covariance matrix
covMatrix <- matrix( c(0.77642931, 0.39590663, 
    0.39590663, 0.49115615), 
    nrow = 2, ncol = 2, byrow = TRUE)

#Create an MxData object including that covariance matrix
data <- mxData(covMatrix, 'cov', numObs = 100)

model <- mxModel(data) 

}
