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

\name{mxEvaluateOnGrid}
\alias{mxEvaluateOnGrid}

\title{Evaluate an algebra on an abscissa grid and collect column results}

\description{
  This function evaluates an algebra on a grid of points provided in
  an auxiliary abscissa matrix.
}

\usage{
mxEvaluateOnGrid(algebra, abscissa)
}


\arguments{
   \item{algebra}{the name of the single column matrix to be evaluated.}
   \item{abscissa}{the name of the abscissa matrix. See details.}
}

\details{
  The abscissa matrix must be in a specific format.
  The variables are in the rows. Abscissa row names must match names of
  free variables.
  The grid points are in columns.
  For each point (column), the free variables are set to the given
  values and the algebra is re-evaluated. The resulting columns are
  collected as the result.
}

\value{
Returns the collected columns.
}

\examples{
library(OpenMx)

test2 <- mxModel("test2",
	mxMatrix(values=1.1, nrow=1, ncol=1, free=TRUE, name="thang"),
	mxMatrix(nrow=1, ncol=1, labels="abscissa1", free=TRUE, name="currentAbscissa"),
	mxMatrix(values=-2:2, nrow=1, ncol=5, name="abscissa",
		 dimnames=list(c('abscissa1'), NULL)),
	mxAlgebra(rbind(currentAbscissa + thang, currentAbscissa * thang), name="stuff"),
	mxAlgebra(mxEvaluateOnGrid(stuff, abscissa), name="grid"))
	
test2 <- mxRun(test2)
omxCheckCloseEnough(test2$grid$result, matrix(c(-1:3 + .1, -2:2 * 1.1), ncol=5, nrow=2,byrow=TRUE))
}
