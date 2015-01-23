%
%   Copyright 2007-2015 The OpenMx Project
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

\name{MxData-class}
\alias{MxData}
\alias{MxNonNullData-class}
\alias{print,MxNonNullData-method}
\alias{show,MxNonNullData-method}

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
    type \tab - \tab Either 'raw', 'cov', or  'cor' \cr
    numObs \tab - \tab The number of oberservations \cr
  }
  
The 'name' slot is the name of the MxData object. %Use of MxData objects in other functions in the \link{OpenMx} library may require reference by name.

The \sQuote{observed} slot is used to contain data, either as a matrix or as a data frame. Use of the data in this slot by other functions depends on the value of the 'type' slot. When 'type' is equal to 'cov' or 'cor', the data input into the 'matrix' slot should be a symmetric matrix or data frame.

The 'vector' slot is used to contain a vector of numeric values, which is used as a vector of means for MxData objects with 'type' equal to 'cov' or 'cor'. This slot may be used in estimation using the \link{mxFitFunctionML} function. 

The 'type' slot may take one of four supported values:
\describe{
\item{raw}{The contents of the \sQuote{observed} slot are treated as raw data. Missing values are permitted and must be designated as the system missing value. The 'vector' and 'numObs' slots cannot be specified, as the 'vector' argument is not relevant and the 'numObs' argument is automatically populated with the number of rows in the data. Data of this type may use the \link{mxFitFunctionML} function as its fit function in MxModel objects, which can deal with covariance estimation under full-information maximum likelihood.}

\item{cov}{The contents of the \sQuote{observed} slot are treated as a covariance matrix. The 'vector' argument is not required, but may be included for estimations involving means. The 'numObs' slot is required. Data of this type may use fit functions such as the \link{mxFitFunctionML}, depending on the specified model.}

\item{cor}{The contents of the \sQuote{observed} slot are treated as a correlation matrix. The 'vector' argument is not required, but may be included for estimations involving means. The 'numObs' slot is required.  Data of this type may use fit functions such as the \link{mxFitFunctionML}, depending on the specified model.}

}

The 'numObs' slot describes the number of observations in the data. If 'type' equals 'raw', then 'numObs' is automatically populated as the number of rows in the matrix or data frame in the \sQuote{observed} slot. If 'type' equals 'cov' or 'cor', then this slot must be input using the 'numObs' argument in the \link{mxData} function when the MxData argument is created.

MxData objects may not be included in \link{MxAlgebra} objects or use the \link{mxFitFunctionAlgebra} function. If these capabilities are desired, data should be appropriately input or transformed using the \link{mxMatrix} and \link{mxAlgebra} functions.

While column names are stored in the \sQuote{observed} slot of MxData objects, these names are not recognized as variable names in \link[=MxPath-class]{MxPath} objects. Variable names must be specified using the 'manifestVars' argument of the \link{mxModel} function prior to use in \link[=MxPath-class]{MxPath} objects.

The mxData function does not currently place restrictions on the size, shape, or symmetry of matrices input into the \sQuote{observed} argument. While it is possible to specify MxData objects as covariance or correlation matrices that do not have the properties commonly associated with these matrices, failure to correctly specify these matrices will likely lead to problems in model estimation.
}

\references{
The OpenMx User's guide can be found at http://openmx.psyc.virginia.edu/documentation.
}

\seealso{
\link{mxData} for creating MxData objects, \link{matrix} and \link{data.frame} for objects which may be entered as arguments in the 'matrix' slot. More information about the OpenMx package may be found \link[=OpenMx]{here}.
}
