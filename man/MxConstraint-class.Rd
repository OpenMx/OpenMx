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

\name{MxConstraint-class}
\alias{MxConstraint-class}
\alias{MxConstraint}

\title{MxConstraint Class}

\description{
   MxConstraint is an S4 class. An MxConstraint object 
   is a \link[=Named-entity]{named entity}. New instances of this class can
   be created using the function \link{mxConstraint}.
}

\details{
The MxConstraint class has the following slots:

  \tabular{rcl}{
    \tab \tab \cr
    name \tab - \tab The name of the object \cr
    formula \tab - \tab The R expression to be evaluated
  }
  
The \sQuote{name} slot is the name of the MxConstraint object. Use of MxConstraint objects in other functions in the \link{OpenMx} library may require reference by name.

The \sQuote{formula} slot is an expression containing the expression to be evaluated. These objects are operated on or related to one another using one or more operations detailed in the \code{\link{mxConstraint}} help file.

Slots may be referenced with the $ symbol. See the documentation for \link[methods]{Classes} and the examples in the \link{mxConstraint} document for more information.
}

\references{
The OpenMx User's guide can be found at http://openmx.psyc.virginia.edu/documentation.
}

\seealso{
\link{mxConstraint} for the function that creates MxConstraint objects.
}
