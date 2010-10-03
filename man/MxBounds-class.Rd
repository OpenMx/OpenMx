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

\name{MxBounds-class}
\alias{MxBounds}

\title{MxBounds Class}

\description{
   MxBounds is an S4 class. New instances of this class can
   be created using the function \link{mxBounds}.
}

\details{
  The MxBounds class has the following slots:

  \tabular{rcl}{
	\tab \tab \cr
	min \tab - \tab The lower bound \cr 
	max \tab - \tab The upper bound \cr
	parameters \tab - \tab The vector of parameter names \cr     
  }
  
The 'min' and 'max' slots hold scalar numeric values for the lower and upper bounds on the list of parameters, respectively.

Parameters may be any free parameter or parameters from an \link{MxMatrix} object. Parameters may be referenced either by name or by referring to their position in the 'spec' matrix of an \code{MxMatrix} object. To affect an estimation or optimization, an MxBounds object must be included in an \link{MxModel} object with all referenced \link{MxAlgebra} and \link{MxMatrix} objects.

Slots may be referenced with the @ symbol. See the documentation for \link[methods]{Classes} and the examples in the \link{mxBounds} document for more information.
}

\references{
The OpenMx User's guide can be found at http://openmx.psyc.virginia.edu/documentation.
}

\seealso{
\link{mxBounds} for the function that creates MxBounds objects. \link{MxMatrix} and \link{mxMatrix} for free parameter specification. More information about the OpenMx package may be found \link[=OpenMx]{here}. 
}
