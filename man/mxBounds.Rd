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

\name{mxBounds}
\alias{mxBounds}

\title{Create MxBounds Object}

\description{
This function creates a new \link{MxBounds} object.
}

\usage{
mxBounds(parameters, min = NA, max = NA)
}

\arguments{
\item{parameters}{A character vector indicating the names of the parameters on which to apply bounds.}
\item{min}{A numeric value for the lower bound. NA means use default value.}
\item{max}{A numeric value for the upper bound. NA means use default value.}
}

\details{
Creates a set of boundaries or limits for a parameter or set of parameters. Parameters may be any free parameter or parameters from an \link{MxMatrix} object. Parameters may be referenced either by name or by referring to their position in the 'spec' matrix of an \code{MxMatrix} object.

Minima and maxima may be specified as scalar numeric values.
}
    
\value{
Returns a new \link{MxBounds} object. If used as an argument in an \link{MxModel} object, the parameters referenced in the 'parameters' argument must also be included prior to optimization.
}

\references{
The OpenMx User's guide can be found at http://openmx.psyc.virginia.edu/documentation.}

\seealso{
\link{MxBounds} for the S4 class created by mxBounds. \link{MxMatrix} and \link{mxMatrix} for free parameter specification. More information about the OpenMx package may be found \link[=OpenMx]{here}. 
}

\examples{
#Create lower and upper bounds for parameters 'A' and 'B'
bounds <- mxBounds(c('A', 'B'), 3, 5)

#Create a lower bound of zero for a set of variance parameters
varianceBounds <- mxBounds(c('Var1', 'Var2', 'Var3'), 0)
}
