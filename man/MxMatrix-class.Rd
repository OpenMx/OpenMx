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

\name{MxMatrix-class}
\alias{MxMatrix}
    
\title{MxMatrix Class}

\description{
   MxMatrix is an S4 class. An MxMatrix object is a \link[=Named-entity]{named entity}.
   New instances of this class can be created using the function \link{mxMatrix}.
   MxMatrix objects may be used as arguments in other functions from the OpenMx library, including \link{mxAlgebra}, \link{mxConstraint}, and \link{mxModel}.
}
    
    \details{
The MxMatrix class has the following slots:
      \tabular{rcl}{
            name \tab - \tab the name of the object \cr
            free \tab - \tab  the free matrix \cr
            values \tab - \tab the values matrix \cr
            labels \tab - \tab  the labels matrix \cr
      }
      
The 'name' slot is the name of the MxMatrix object. Use of MxMatrix objects in an \link{mxAlgebra} or \link{mxConstraint} function requires reference by name.

The 'free' slot takes a matrix which describes the location of free and fixed parameters. A variable is a free parameter if-and-only-if the corresponding value in the 'free' matrix is 'TRUE'. Free parameters are elements of an MxMatrix object whose values may be changed by an objective function when that MxMatrix object is included in an \link{MxModel} object and evaluated using the \link{mxRun} function.

The 'values' slot takes a matrix of numeric values. If an element is specified as a fixed parameter in the 'free' matrix, then the element in the 'values' matrix is treated as a constant value and cannot be altered or updated by an objective function when included in an \link{mxRun} function. If an element is specified as a free parameter in the 'free' matrix, the element in the 'value' matrix is considered a starting value and can be changed by an objective function when included in an \link{mxRun} function. 

The 'labels' slot takes a matrix which describes the labels of free and fixed parameters. Fixed parameters with identical labels must have identical values.  Free parameters with identical labels impose an equality constraint. The same label cannot be applied to a free parameter and a fixed parameter. A free parameter with the label 'NA' implies a unique free parameter, that cannot be constrained to equal any other free parameter.

}

\references{
The OpenMx User's guide can be found at http://openmx.psyc.virginia.edu/documentation.
}

\seealso{
\link{mxMatrix} for creating MxMatrix objects. More information about the OpenMx package may be found \link[=OpenMx]{here}. 
}
