%
%   Copyright 2007-2012 The OpenMx Project
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

\name{omxMatrixOperations}
\alias{omxMatrixOperations}
\alias{omxCbind}
\alias{omxRbind}
\alias{omxTranspose}

\title{MxMatrix operations}

\description{
   \code{omxCbind} columnwise binding of two or more MxMatrices.
   \code{omxRbind} rowwise binding of two or more MxMatrices.
   \code{omxTranspose} transpose of MxMatrix.
}

\usage{

omxCbind(..., allowUnlabeled = 
	getOption("mxOptions")[["Allow Unlabeled"]], 
	dimnames = NA, name = NA)
omxRbind(..., allowUnlabeled = 
	getOption("mxOptions")[["Allow Unlabeled"]], 
	dimnames = NA, name = NA)
omxTranspose(matrix, allowUnlabeled = 
	getOption("mxOptions")[["Allow Unlabeled"]], 
	dimnames = NA, name = NA)
}

\arguments{
   \item{...}{two or more MxMatrix objects}
   \item{allowUnlabeled}{whether or not to accept free parameters with NA labels}
   \item{dimnames}{list. The dimnames attribute for the matrix: a list of length 2 giving the row and column names respectively. An empty list is treated as NULL, and a list of length one as row names. The list can be named, and the list names will be used as names for the dimensions.}
   \item{name}{an optional character string indicating the name of the MxMatrix object}
}
