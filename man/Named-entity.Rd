%
%   Copyright 2007-2014 The OpenMx Project
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

\name{Named-entity}
\alias{Named-entity}
\alias{Named-entities}
\alias{named-entity}
\alias{named-entities}

\title{Named Entities}

\description{
    A named entity is an S4 object that can be referenced by name.
}

\details{
Every named entity is guaranteed to have a slot called "name". Within a model, the named entities of that model can be accessed using the $ operator.  Access is limited to one nesting depth, such that if 'B' is a submodel of 'A', and 'C' is a matrix of 'B', then 'C' must be accessed using A$B$C. 

The following S4 classes are named entities in the OpenMx library: \link{MxAlgebra}, \link{MxConstraint},  \link{MxMatrix}, \link{MxModel}, \link{MxData}, and MxObjective.
}

\examples{

library(OpenMx)

# Create a model, add a matrix to it, and then access the matrix by name.

testModel <- mxModel(model="anEmptyModel")

testMatrix <- mxMatrix(type="Full", nrow=2, ncol=2, values=c(1,2,3,4), name="yourMatrix")

yourModel <- mxModel(testModel, testMatrix, name="noLongerEmpty")

yourModel$yourMatrix



}
