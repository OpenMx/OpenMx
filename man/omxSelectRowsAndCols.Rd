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

\name{omxSelectRowsAndCols}
\alias{omxSelectRowsAndCols}
\alias{omxSelectRows}
\alias{omxSelectCols}

\title{Filter rows and columns from an mxMatrix}

\description{
   This function filters rows and columns from a matrix using a single row or column R matrix as a selector. 
}

\usage{
omxSelectRowsAndCols(x, selector)
omxSelectRows(x, selector)
omxSelectCols(x, selector)
}


\arguments{
   \item{x}{the matrix to be filtered}
   \item{selector}{A single row or single column R matrix indicating which values should be filtered from the mxMatrix.}
}

\details{

omxSelectRowsAndCols, omxSelectRows, and omxSelectCols returns the filtered entries in a target matrix specified by a single row or single column selector matrix. Each entry in the selector matrix is treated as a logical data indicating if the corresponding entry in the target matrix should be excluded (0 or FALSE) or included (not 0 or TRUE). Typically the function is used to filter data from a target matrix using an existence vector which specifies what data entries are missing. This can be seen in the demo: RowObjectiveFIMLBivariateSaturated. 
}

\value{
Returns a new matrix with the filtered data.
}

\references{
The function is most often used when filtering data for missingness. This can be seen in the demo: RowObjectiveFIMLBivariateSaturated. The OpenMx User's guide can be found at http://openmx.psyc.virginia.edu/documentation. The omxSelect* functions share some similarity to the Extract function in the R programming language.
}

\examples{
loadings <- c(1, -0.625, 0.1953125, 1, -0.375, 0.0703125, 1, -0.375, 0.0703125)
loadings <- matrix(loadings, 3, 3, byrow= TRUE)
existenceList <- c(1, 0, 1)
existenceList <- matrix(existenceList, 1, 3, byrow= TRUE)
rowsAndCols <- omxSelectRowsAndCols(loadings, existenceList)
rows <- omxSelectRows(loadings, existenceList)
cols <- omxSelectCols(loadings, existenceList)
}
