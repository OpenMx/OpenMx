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

\name{mxExpectationStateSpace}
\alias{mxExpectationStateSpace}

\title{Create an MxExpectationStateSpace Object}

\description{
This function creates a new MxExpectationStateSpace object.
}

\usage{
mxExpectationStateSpace(A="A", B="B", C="C", D="D", Q="Q", R="R", x0="x0", P0="P0", 
                        dimnames = NA, thresholds = NA, threshnames = dimnames)
}


\arguments{
   \item{A}{A character string indicating the name of the 'A' matrix.}
   \item{B}{A character string indicating the name of the 'B' matrix.}
   \item{C}{A character string indicating the name of the 'C' matrix.}
   \item{D}{A character string indicating the name of the 'D' matrix.}
   \item{Q}{A character string indicating the name of the 'Q' matrix.}
   \item{R}{A character string indicating the name of the 'R' matrix.}
   \item{x0}{A character string indicating the name of the 'x0' matrix.}
   \item{P0}{A character string indicating the name of the 'P0' matrix.}
   \item{dimnames}{An optional character vector to be assigned to the column names of the 'F' and 'M' matrices.}   
   \item{thresholds}{An optional character string indicating the name of the thresholds matrix.}
   \item{threshnames}{An optional character vector to be assigned to the column names of the thresholds matrix.}
}

\details{
Expectation functions define the way that model expectations are calculated. The mxExpectationStateSpace provides maximum likelihood estimates of free parameters in a model of the covariance of a given \link{MxData} object. This model is defined by reticular action modeling (McArdle and McDonald, 1984). The 'A', 'B', 'C', 'D', 'Q', 'R', 'x0', and 'P0' arguments must refer to \link{MxMatrix} objects with the associated properties of the A, B, C, D, Q, R, x0, and P0 matrices in the RAM modeling approach.

The state space expectation is defined by the following equations.

\deqn{x_{t+1} = A x_t + B u_t + q_t}
\deqn{y_t = C x_t + D u_t + r_t}

The 'dimnames' arguments takes an optional character vector.  If this argument is not a single NA, then this vector be assigned to be the column names of the '?' matrix and optionally to the '?' matrix, if the '?' matrix exists.

THE NEXT PART IS INCOMPLETE.  MIKE HUNTER NEEDS TO FINISH THIS PART.

The 'A' argument refers to the A matrix in the State Space approach. This matrix consists of . 

The 'B' argument refers to the A matrix in the State Space approach. This matrix consists of . 

The 'C' argument refers to the A matrix in the State Space approach. This matrix consists of . 

The 'D' argument refers to the A matrix in the State Space approach. This matrix consists of . 

The 'Q' argument refers to the A matrix in the State Space approach. This matrix consists of . 

The 'R' argument refers to the A matrix in the State Space approach. This matrix consists of . 

The 'x0' argument refers to the A matrix in the State Space approach. This matrix consists of . 

The 'P0' argument refers to the A matrix in the State Space approach. This matrix consists of . 


The \link{MxMatrix} objects included as arguments may be of any type, but should have the properties described above. The mxExpectationStateSpace will not return an error for incorrect specification, but incorrect specification will likely lead to estimation problems or errors in the \link{mxRun} function.

mxExpectationStateSpace evaluates with respect to an \link{MxData} object. The \link{MxData} object need not be referenced in the mxExpectationStateSpace function, but must be included in the \link{MxModel} object. mxExpectationStateSpace requires that the 'type' argument in the associated \link{MxData} object be equal to 'cov' or 'cor'.

To evaluate, place mxExpectationStateSpace objects, the \link{mxData} object for which the expected covariance approximates, referenced \link{MxAlgebra} and \link{MxMatrix} objects, and optional \link{MxBounds} and \link{MxConstraint} objects in an \link{MxModel} object. This model may then be evaluated using the \link{mxRun} function. The results of the optimization can be found in the 'output' slot of the resulting model, and may be obtained using the \link{mxEval} function..
}

\value{
    Returns a new MxExpectationStateSpace object. mxExpectationStateSpace objects should be included with models with referenced \link{MxAlgebra}, \link{MxData} and \link{MxMatrix} objects.
}

\references{
McArdle, J. J. and MacDonald, R. P. (1984). Some algebraic properties of the Reticular Action Model for moment structures. \emph{British Journal of Mathematical and Statistical Psychology, 37,} 234-251.

The OpenMx User's guide can be found at http://openmx.psyc.virginia.edu/documentation.
}

\examples{
   
# Create and fit a model using mxMatrix, mxAlgebra, mxExpectationNormal, and mxFitFunctionML

library(OpenMx)

# Simulate some data

# The example is that there is no example.

}
