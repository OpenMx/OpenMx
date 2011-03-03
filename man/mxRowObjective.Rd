%
%   Copyright 2007-2010 The OpenMx Project
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

% Author: Michael D. Hunter
% Revision History
%   Fri Oct 8 15:56:13 EDT 2010 - Created file
%                                 Modified from mxFIMLObjective.Rd


\name{mxRowObjective}
\alias{mxRowObjective}

\title{Create MxRowObjective Object}

\description{
   This function creates a new MxRowObjective object.
}

\usage{
mxRowObjective(rowAlgebra, reduceAlgebra, dimnames,
   rowResults = "rowResults", filteredDataRow = "filteredDataRow",
   existenceVector = "existenceVector")
}


\arguments{
   \item{rowAlgebra}{A character string indicating the name of the algebra to be evaluated row-wise.}
   \item{reduceAlgebra}{A character string indicating the name of the algebra that collapses the row results into a single number which is then optimized.}
   \item{dimnames}{A character vector of names corresponding to columns be extracted from the data set.}
   \item{rowResults}{The name of the auto-generated "rowResults" matrix.  See details.}
   \item{filteredDataRow}{The name of the auto-generated "filteredDataRow" matrix. See details.}
   \item{existenceVector}{The name of the auto-generated "existenceVector" matrix. See details.}
}

\details{
Objective functions are functions for which free parameter values are chosen such that the value of the objective function is minimized. The mxRowObjective function evaluates a user-defined \link{MxAlgebra} object called the \sQuote{rowAlgebra} in a row-wise fashion.  It then stores results of the row-wise evaluation in another \link{MxAlgebra} object called the \sQuote{rowResults}.  Finally, the mxRowObjective function collapses the row results into a single number which is then used for optimization.  The \link{MxAlgebra} object named by the \sQuote{reduceAlgebra} collapses the row results into a single number.

The \sQuote{filteredDataRow} is populated in a row-by-row fashion with all the non-missing data from the current row.  You cannot assume that the length of the filteredDataRow matrix remains constant (unless you have no missing data). The \sQuote{existenceVector} is populated in a row-by-row fashion with a value of 1.0 in column j if a non-missing value is present in the data set in column j, and a value of 0.0 otherwise. Use the functions \link{omxSelectRowRows}, \link{omxSelectCols}, and \link{omxSelectRowsAndCols} to shrink other matrices so that their dimensions will be conformable to the size of \sQuote{filteredDataRow}. 
}

\value{
Returns a new MxRowObjective object. MxRowObjective objects should be included with models with referenced \link{MxAlgebra} objects.
}

\references{
The OpenMx User's guide can be found at http://openmx.psyc.virginia.edu/documentation.
}

\examples{

# Model that adds two data columns row-wise, then sums that column
# Notice no optimization is performed here.
xdat <- data.frame(a=rnorm(10), b=1:10) # Make data set
amod <- mxModel(
		    mxData(observed=xdat, type='raw'),
		    mxAlgebra(sum(filteredDataRow), name = 'rowAlgebra'),
		    mxAlgebra(sum(rowResults), name = 'reduceAlgebra'),
		    mxRowObjective(
			    rowAlgebra='rowAlgebra',
			    reduceAlgebra='reduceAlgebra',
			    dimnames=c('a','b'))
)

# Model that find the parameter that minimizes the sum of the
#  squared difference between the parameter and a data row.
bmod <- mxModel(
		    mxData(observed=xdat, type='raw'),
		    mxMatrix(values=.75, ncol=1, nrow=1, free=TRUE, name='B'),
		    mxAlgebra((filteredDataRow - B) ^ 2, name='rowAlgebra'),
		    mxAlgebra(sum(rowResults), name='reduceAlgebra'),
		    mxRowObjective(
			    rowAlgebra='rowAlgebra',
			    reduceAlgebra='reduceAlgebra',
			    dimnames=c('a'))
)


}
