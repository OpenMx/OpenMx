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

\name{MxCI-class}
\alias{MxCI}

\title{MxCI Class}

\description{
   MxCI is an S4 class. An MxCI object is a \link[=Named-entity]{named entity}.
   New instances of this class can be created using the function \link{mxCI}.
   MxCI objects may be used as arguments in the \link{mxModel} function.
}

\details{
The MxCI class has the following slots:

  \tabular{rcl}{
    \tab \tab \cr
    reference \tab - \tab The name of the object \cr
    lowerdelta \tab - \tab Either a matrix or a data frame \cr
    upperdelta \tab - \tab A vector for means, or NA if missing \cr
  }
  
The reference slot contains a character vector of named free parameters, \link[=MxMatrix]{MxMatrices} and \link[=MxAlgebra]{MxAlgebras} on which confidence intervals are desired. Individual elements of \link[=MxMatrix]{MxMatrices} and \link[=MxAlgebra]{MxAlgebras} may be listed as well, using the syntax \dQuote{matrix[row,col]} (see \link{Extract} for more information).

The lowerdelta and upperdelta slots give the changes in likelihoods used to define the confidence interval. The upper bound of the likelihood-based confidence interval is estimated by increasing the parameter estimate, leaving all other parameters free, until the model -2 log likelihood increased by \sQuote{upperdelta}. The lower bound of the confidence interval is estimated by decreasing the parameter estimate, leaving all other parameters free, until the model -2 log likelihood increased by \sQuote{lowerdata}.

Likelihood-based confidence intervals may be specified by including one or more MxCI objects in an \link{MxModel} object. Estimation of confidence intervals requires model optimization using the \link{mxRun} function with the \sQuote{intervals} argument set to TRUE. The calculation of likelihood-based confidence intervals can be computationally intensive, and may add a significant amount of time to model estimation when many confidence intervals are requested.

}

\references{
The OpenMx User's guide can be found at http://openmx.psyc.virginia.edu/documentation.
}

\seealso{
\link{mxCI} for creating MxCI objects. More information about the OpenMx package may be found \link[=OpenMx]{here}.
}
