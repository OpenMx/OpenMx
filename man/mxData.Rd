%
%   Copyright 2007-2019 by the individuals mentioned in the source code history
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
   This function creates a new \link{MxData} object. This can be used all forms of analysis (including WLS: see \link{mxFitFunctionWLS}).
   It packages observed data (e.g. a dataframe, matrix, or cov or cor matrix) into an object with additional information allowing it to be processed in an mxModel.
}

\usage{
   mxData(observed, type, means = NA, numObs = NA, acov=NA, fullWeight=NA,
          thresholds=NA, ..., observedStats=NA, sort=NA, primaryKey = as.character(NA),
          weight = as.character(NA), frequency = as.character(NA),
          verbose = 0L, .parallel=TRUE, .noExoOptimize=TRUE,
     minVariance=sqrt(.Machine$double.eps), algebra=c(),
   warnNPDacov=TRUE, exoFree=NULL,
   naAction=c("pass","fail","omit","exclude"))
}

\arguments{
   \item{observed}{A matrix or data.frame which provides data to the MxData object.}
   \item{type}{A character string defining the type of data in the \sQuote{observed} argument. Must be one of \dQuote{raw}, \dQuote{cov}, \dQuote{cor}, or \dQuote{acov}.}
   \item{means}{An optional vector of means for use when \sQuote{type} is \dQuote{cov}, or \dQuote{cor}.}
   \item{numObs}{The number of observations in the data supplied in the \sQuote{observed} argument. Required unless \sQuote{type} equals \dQuote{raw}.}
   \item{...}{Not used. Forces remaining arguments to be specified by name.}
   \item{observedStats}{A list containing observed statistics for weighted least squares estimation. See details for contents}
   \item{sort}{Whether to sort raw data prior to use (default NA).}
   \item{primaryKey}{The column name of the primary key used to uniquely identify rows (default NA)}
   \item{weight}{The column name containing row weights.}
   \item{frequency}{The column name containing row frequencies.}
   \item{verbose}{level of diagnostic output.}
   \item{.parallel}{logical. Whether to compute observed summary statistics in parallel.}
   \item{.noExoOptimize}{logical. Whether to use math short-cuts for the case of no exogenous predictors.}
   \item{minVariance}{numeric. The minimum acceptable variance for \sQuote{observedStats$cov}.}
   \item{acov}{Deprecated in favor of the acov element of observedStats.}
   \item{fullWeight}{Deprecated in favor of the fullWeight element of observedStats.}
   \item{thresholds}{Deprecated in favor of the thresholds element of observedStats.}
   \item{algebra}{character vector. Names of algebras used to fill in
     calculated columns of raw data. \lifecycle{experimental}}
   \item{warnNPDacov}{logical. Whether to warn when the asymptotic
     covariance matrix is non-positive definite.}
   \item{exoFree}{logical matrix of observed manifests by
   exogenous predictors. Defaults to all TRUE, but you can fix some
   regression coefficients in the \code{observedStats} \code{slope}
   matrix to zero by setting entries to FALSE. \lifecycle{experimental}}
 \item{naAction}{Specify treatment of missing data. See details. \lifecycle{maturing}}
}

\details{
The mxData function creates \link{MxData} objects used in \link{mxModel}s.
The \sQuote{observed} argument may take either a data frame or a matrix, which is then described with the
\sQuote{type} argument. Data types describe compatibility and usage with expectation functions in MxModel
objects. Three data types are supported (acov is deprecated).

\describe{
\item{raw}{The contents of the \sQuote{observed} argument are treated as raw data. Missing values are permitted and must be
designated as the system missing value. The \sQuote{means} and \sQuote{numObs} arguments cannot be specified, as the
\sQuote{means} argument is not relevant and the \sQuote{numObs} argument is automatically populated with the number of rows
in the data. Data of this type may use fit functions such as \link{mxFitFunctionML} or \link{mxFitFunctionWLS}.
\link{mxFitFunctionML} will automatically use use full-information maximum likelihood for raw data.}


\item{cov}{The contents of the \sQuote{observed} argument are treated as a covariance matrix. The \sQuote{means} argument is
not required, but may be included for estimations involving means. The \sQuote{numObs} argument is required, which should
reflect the number of observations or rows in the data described by the covariance matrix. Cov data typically use the
\link{mxFitFunctionML} fit function, depending on the specified model.}

\item{acov}{ This type was used for WLS data as created by \link{mxDataWLS}. Unless you are using summary data, its use is deprecated.
Instead, use type =\sQuote{raw} and an \link{mxFitFunctionWLS}. If type \sQuote{acov} is set, the \sQuote{observed} argument will
(usually) contain raw data and the \sQuote{observedStats} slot contain a list of observed statistics.}

\item{cor}{The contents of the \sQuote{observed} argument are treated as a correlation matrix. The \sQuote{means} argument is
not required, but may be included for estimations involving means. The \sQuote{numObs} argument is required, which should
reflect the number of observations or rows in the data described by the covariance matrix. Models with cor data typically use
the \link{mxFitFunctionML} fit function.}
}

\emph{Note on data handling}: OpenMx uses the names of variables to map them onto other elements of your model, such as expectation functions.
Thus for data provided as a data.frame, ensure the columns have appropriate \code{\link{names}}.
Covariance and correlation matrices need to have both the row and column names set and these must be
identical, for instance by using \code{dimnames = list(varNames, varNames)}.

\strong{Correlation data}

To obtain accurate parameter estimates and standard errors,
it is necessary to constrain the model implied covariance matrix to have
unit variances. This constraint is added automatically if you use an
\code{\link{mxModel}} with \code{type='RAM'} or \code{type='LISREL'}.
Otherwise, you will need to add this constraint yourself.

\strong{WLS data}

The \code{observedStats} contains the following named objects: cov, slope, means, acov, fullWeight, and thresholds.

\sQuote{cov} The (polychoric) covariance matrix of raw data variables. An error is raised if any variance is smaller \code{minVariance}.

\sQuote{slope} The regression coefficients from all exogenous predictors to all observed variables. Required for exogenous predictors.

\sQuote{means} The means of the data variables. Required for estimations involving means.

\sQuote{thresholds} Thresholds of ordinal variables. Required for models including ordinal variables.

\sQuote{acov} The asymptotic covariance matrix. Can be diagonal for diagonally weighted least squares.

\sQuote{fullWeight} (optional) is the complete asymptotic weight matrix (all entries non-zero).

\emph{note}: WLS data typically use the \link{mxFitFunctionWLS} function.

\emph{IMPORTANT}: The WLS interface is under heavy development to support both very fast backend processing of raw data while
continuing to support modeling applications which require direct access to the object in the front end. Some user-interface
changes should be expected as we optimize both these workflows.


\strong{Missing values}

For raw data, the \sQuote{naAction} option controls the treatment of missing values.
When set to \sQuote{pass}, the data is passed as-is.
When set to \sQuote{fail}, the presence of any missing value will
trigger an error.
When set to \sQuote{omit}, missing data will be discarded row-wise.
For example, a single missing value in a row will cause the whole row to
be discarded.
When set to \sQuote{exclude}, rows with missing data are retained
but their \sQuote{frequency} is set to zero.

\strong{Weights}

In the case of raw data, the optional \sQuote{weight} argument names a column in the data that contains per-row weights.
Similarly, the optional \sQuote{frequency} argument names a column in the \sQuote{observed} data that contains per-row
frequencies. Frequencies must be integers but weights can be arbitrary real numbers. For data with many repeated response
patterns, organizing the data into unique patterns and frequencies can reduce model evaluation time.

In some cases, the fit function can be evaluated more efficiently when data are sorted. When a primary key is provided,
sorting is disabled. Otherwise, sort defaults to TRUE.

The mxData function does not currently place restrictions on the size, shape, or symmetry of matrices input into the
\sQuote{observed} argument. While it is possible to specify MxData objects as covariance or correlation matrices that do not
have the properties commonly associated with these matrices, failure to correctly specify these matrices will likely lead to
problems in model estimation.

\emph{note}: MxData objects may not be included in \link{mxAlgebra}s nor in the \link{mxFitFunctionAlgebra} function. To
reference data in these functions, use a \link{mxMatrix} or a definition variable (data.var) label.

Also, while column names are stored in the \sQuote{observed} slot of MxData objects, these names are not automatically
recognized as variable names in \link[=MxPath-class]{mxPaths} in RAM models. These models use the \sQuote{manifestVars} of
the \link{mxModel} function to explicitly identify used variables used in the model.
}

\value{
    Returns a new \link{MxData} object.
}

\references{
The OpenMx User's guide can be found at http://openmx.ssri.psu.edu/documentation.
}

\seealso{
To generate data, see \code{\link{mxGenerateData}}; For objects which may be entered as arguments in the
\sQuote{observed} slot, see \link{matrix} and \link{data.frame}. See \link{MxData} for the S4 class created by mxData.
For WLS data, see \link{mxDataWLS} (deprecated). More information about the OpenMx package may be found
\link[=OpenMx]{here}.

 
}

\examples{  
    
library(OpenMx)

# Simple covariance model. See other mxFitFunctions for examples with different data types

# 1. Create a covariance matrix x and y
covMatrix <- matrix(nrow = 2, ncol = 2, byrow = TRUE,
	c(0.77642931, 0.39590663, 
      0.39590663, 0.49115615)
)
covNames <- c("x", "y")
dimList <- list(covNames, covNames)
dimnames(covMatrix) <- dimList

# 2. Create an MxData object from covMatrix
testData <- mxData(observed=covMatrix, type="cov", numObs = 100)

testModel <- mxModel(model="testModel2",
	mxMatrix(name="expCov", type="Symm", nrow=2, ncol=2,
                 values=c(.2,.1,.2), free=TRUE, dimnames=dimList),
    mxExpectationNormal("expCov", dimnames=covNames),
    mxFitFunctionML(),
	testData
)

outModel <- mxRun(testModel)

summary(outModel)

}
