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
mxRowObjective(rowAlgebra, rowResults = NA, reduceAlgebra = NA)
}


\arguments{
   \item{rowAlgebra}{A character string indicating the name of the algebra to be evaluated row-wise.}
   \item{rowResults}{A character string indicating the name of the algebra that stores the results of the row algebra.}
   \item{reduceAlgebra}{A character string indicating the name of the algebra that collapses the row results into a single number which is then optimized.}
}

\details{
Objective functions are functions for which free parameter values are chosen such that the value of the objective function is minimized. The mxRowObjective function evaluates a user-defined \link{MxAlgebra} object called the 'rowAlgebra' in a row-wise fashion.  It then stores results of the row-wise evaluation in another \link{MxAlgebra} object called the 'rowResults'.  Finally, the mxRowObjective function collapses the row results into a single number which is then used for optimization.  The \link{MxAlgebra} object named by the 'reduceAlgebra' collapses the row results into a single number.  Any free parameters found in the 'rowAlgebra' will be optimized to minimize the result of the 'reduceAlgebra'.

All arguments to the mxRowObjective function are required.  There are no optional arguments.

mxRowObjective evaluates with respect to several \link{MxAlgebra} objects. The \link{MxAlgebra} objects need to be referenced in the mxRowObjective function, and also must be included in the \link{MxModel} object.  Missing values are not currently permitted in the associated \link{MxAlgebra} objects.  \link{MxData} objects can be indirectly referenced by mxRowObjective by placing all or part of the \link{MxData} object into an \link{MxAlgebra} object and referencing the \link{MxAlgebra} object.

To evaluate, place MxRowObjective objects, referenced \link{MxAlgebra} objects, and optional \link{MxData}, \link{MxMatrix}, \link{MxBounds}, and \link{MxConstraint} objects in an \link{MxModel} object. This model may then be evaluated using the \link{mxRun} function. The results of the optimization can be found in the 'output' slot of the resulting model, and may be referenced using the \link{Extract} functionality.
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
