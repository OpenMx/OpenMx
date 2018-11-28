%
%   Copyright 2007-2018 by the individuals mentioned in the source code history
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
   This function creates a new \link{MxData} object. (for WLS data see \link{mxDataWLS}).
}

\usage{
   mxData(observed, type, means = NA, numObs = NA, acov=NA, fullWeight=NA,
          thresholds=NA, ..., observedStats=NA, sort=NA, primaryKey = as.character(NA),
          weight = as.character(NA), frequency = as.character(NA),
          verbose = 0L, preferredFit = "ML")
}

\arguments{
   \item{observed}{A matrix or data.frame which provides data to the MxData object.}
   \item{type}{A character string defining the type of data in the \sQuote{observed} argument. Must be one of \dQuote{raw}, \dQuote{cov}, \dQuote{cor}, \dQuote{'sscp'} or \dQuote{acov}.}
   \item{means}{An optional vector of means for use when \sQuote{type} is \dQuote{cov}, or \dQuote{cor}.}
   \item{numObs}{The number of observations in the data supplied in the \sQuote{observed} argument. Required unless \sQuote{type} equals \dQuote{raw}.}
   \item{acov}{Deprecated.}
   \item{fullWeight}{Deprecated.}
   \item{thresholds}{Deprecated.}
   \item{...}{Not used. Forces remaining arguments to be specified by name.}
   \item{observedStats}{A list containing observed statistics for
   weighted least squares estimation.}
   \item{sort}{Whether to sort raw data prior to use (default NA)}
   \item{primaryKey}{The column name of the primary key used to uniquely
     identify rows (default NA)}
   \item{weight}{The column name containing row weights}
   \item{frequency}{The column name containing row frequencies}
\item{verbose}{level of diagnostic output}
\item{preferredFit}{used to keep track of whether the MxData object was
   created by mxDataWLS}
}

\details{
The mxData function creates \link{MxData} objects, which can be used as arguments in \link{MxModel} objects. The \sQuote{observed} argument may take either a data frame or a matrix, which is then described with the \sQuote{type} argument. Data types describe compatibility and usage with expectation functions in MxModel objects. Four different data types are supported (a fifth, sscp, is not yet implemented):

\describe{
\item{raw}{The contents of the \sQuote{observed} argument are treated as raw data. Missing values are permitted and must be designated as the system missing value. The \sQuote{means} and \sQuote{numObs} arguments cannot be specified, as the \sQuote{means} argument is not relevant and the \sQuote{numObs} argument is automatically populated with the number of rows in the data. Data of this type may use fit functions such as \link{mxFitFunctionML} function in MxModel objects, which will automatically use covariance estimation under full-information maximum likelihood for this data type.}

\item{cov}{The contents of the \sQuote{observed} argument are treated as a covariance matrix. The \sQuote{means} argument is not required, but may be included for estimations involving means. The \sQuote{numObs} argument is required, which should reflect the number of observations or rows in the data described by the covariance matrix. Data of this type may use the fit functions such as \link{mxFitFunctionML}, depending on the specified model.}

\item{cor}{The contents of the \sQuote{observed} argument are treated as a correlation matrix. The \sQuote{means} argument is not required, but may be included for estimations involving means. The \sQuote{numObs} argument is required, which should reflect the number of observations or rows in the data described by the covariance matrix. Data of this type may use the fit functions such as \link{mxFitFunctionML} functions, depending on the specified model.}

\item{acov}{The best way to create data of this type is to use the
  \link{mxDataWLS} function.  The contents of the \sQuote{observed}
  argument is regarded as raw data. The \sQuote{observedStats} slot
  contains a list of observed statistics. The list is expected to
  contain a number of named objects: \sQuote{cov} is the (polychoric)
  covariance matrix of the continuous and ordinal
  variables. \sQuote{means} (optional) may be included for estimations
  involving means.  \sQuote{thresholds} should be included for
  estimation involving thresholds and ordinal variables. \sQuote{acov}
  is the asymptotic covariance matrix. \sQuote{acov} can be diagonal for
  daigonally weighted least squares. \sQuote{fullWeight} (optional) is
  the complete aymptotic weight matrix (all entries non-zero). Data of this
  type almost certainly use the \link{mxFitFunctionWLS} function, but
  may depend on the specified model.}
}

MxData objects may not be included in \link{MxAlgebra} objects or use the \link{mxFitFunctionAlgebra} function. If these capabilities are desired, data should be appropriately input or transformed using the \link{mxMatrix} and \link{mxAlgebra} functions.

While column names are stored in the \sQuote{observed} slot of MxData objects, these names are not recognized as variable names in \link[=MxPath-class]{MxPath} objects. Variable names must be specified using the \sQuote{manifestVars} argument of the \link{mxModel} function prior to use in \link[=MxPath-class]{MxPath} objects.

The mxData function does not currently place restrictions on the size, shape, or symmetry of matrices input into the \sQuote{observed} argument. While it is possible to specify MxData objects as covariance or correlation  matrices that do not have the properties commonly associated with these matrices, failure to correctly specify these matrices will likely lead to problems in model estimation.

OpenMx uses the names of variables to map them onto the expectation functions and other elements associated with your model. For data.frames, ensure you have set the names(). For matrices set names using, for instance, row.names=c(\dQuote{your}, \dQuote{columns}). Covariance and correlation matrices need to have both the row and column names set and these must be identical, for instance by using \code{dimnames = list(varNames, varNames)}.

In the case of raw data, the optional \sQuote{weight} argument names a
column in the \sQuote{observed} data that contains per-row weights.
Similarly, the optional \sQuote{frequency} argument names a column
in the \sQuote{observed} data that contains per-row frequencies.
Frequencies must be integral but weights can be arbitrary real numbers.
For data with many repeated response patterns, organizing the data into
unique patterns and frequencies can reduce model evaluation time.

In some cases, the fit function can be evaluated more efficiently when
data are sorted.
When a primary key is provided, sorting is disabled. Otherwise, sort defaults to TRUE.
}

\value{
    Returns a new \link{MxData} object.
}

\references{
The OpenMx User's guide can be found at http://openmx.ssri.psu.edu/documentation.
}

\seealso{
For WLS data, see \link{mxDataWLS}. To generate data, see \code{\link{mxGenerateData}}; For objects which may be entered as arguments in the \sQuote{observed} slot, see \link{matrix} and \link{data.frame}. See \link{MxData} for the S4 class created by mxData. More information about the OpenMx package may be found \link[=OpenMx]{here}. 
}

\examples{  
    
library(OpenMx)

# Create a covariance matrix
covMatrix <- matrix( c(0.77642931, 0.39590663, 
    0.39590663, 0.49115615), 
    nrow = 2, ncol = 2, byrow = TRUE)
covNames <- c("x", "y")
dimList <- list(covNames, covNames)
dimnames(covMatrix) <- dimList

# Create an MxData object including that covariance matrix
testData <- mxData(observed=covMatrix, type="cov", numObs = 100)

testModel <- mxModel(model="testModel",
                mxMatrix(type="Symm", nrow=2, ncol=2, values=c(.2,.1,.2), 
                         free=TRUE, name="expCov", dimnames=dimList),
                mxExpectationNormal(covariance="expCov", dimnames=covNames),
                mxFitFunctionML(),
                testData) 

outModel <- mxRun(testModel)

summary(outModel)

}
