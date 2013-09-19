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
Expectation functions define the way that model expectations are calculated. When used in conjunction with the \link{mxFitFunctionML}, the mxExpectationStateSpace uses maximum likelihood prediction error decomposition (PED) to obtain estimates of free parameters in a model of the raw \link{MxData} object. State space expectations treat the raw data as a multivariate time series of equally spaced times with each row corresponding to a single occasion.  This is not a model of the block Toeplitz lagged autocovariance matrix.  State space expectations implement a classical Kalman filter to produce expectations.

The following are not yet implemented: square root Kalman filter (in Cholesky or singular value decomposition form), extended Kalman filter for linear approximations to nonlinear state space models, unscented Kalman filter for highly nonlinear state space models, Kalman-Bucy filter for continuous time modeling, hybrid Kalman filter for continuous latent time with discrete observations, and Rauch-Tung-Striebel smoother for updating forecast state estimates after a complete forward pass through the data has been made.

This model uses notation for the model matrices commonly found in engineering and control theory.  

 The 'A', 'B', 'C', 'D', 'Q', 'R', 'x0', and 'P0' arguments must be the names of \link{MxMatrix} objects with the associated properties of the A, B, C, D, Q, R, x0, and P0 matrices in the state space modeling approach.

The state space expectation is defined by the following model equations.

\deqn{x_{t+1} = A x_t + B u_t + q_t}
\deqn{y_t = C x_t + D u_t + r_t}

The first equation is called the state equation.  It describes how the latent states change over time.  Also, the state equation in state space modeling is directly analogous to the structural model in LISREL structural equation modeling.

The second equation is called the output equation.  It describes how the latent states relate to the observed states at a single point in time.  The output equation shows how the observed output is produced by the latent states.  Also, the output equation in state space modeling is directly analogous to the measurement model in LISREL structural equation modeling.

The state and output equations, together with some minimal assumptions, imply ... covariance and means

The 'dimnames' arguments takes an optional character vector.  % TODO: Decide how dimnames are handled.
%If this argument is not a single NA, then this vector be assigned to be the row names of the 'C' matrix and optionally to the '?' matrix, if the '?' matrix exists.

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
   
# Create and fit a model using mxMatrix, mxExpectationStateSpace, and mxFitFunctionML

require(OpenMx)
data(demoOneFactor)
nvar <- ncol(demoOneFactor)
ssModel <- mxModel(name="State Space Manual Example",
	mxMatrix("Full", 1, 1, TRUE, .3, name="A"),
	mxMatrix("Zero", 1, 1, name="B"),
	mxMatrix("Full", nvar, 1, TRUE, .6, name="C", dimnames=list(colnames(demoOneFactor)
, "F1")),
	mxMatrix("Zero", nvar, 1, name="D"),
	mxMatrix("Diag", 1, 1, FALSE, 1, name="Q"),
	mxMatrix("Diag", nvar, nvar, TRUE, .2, name="R"),
	mxMatrix("Zero", 1, 1, name="x0"),
	mxMatrix("Diag", 1, 1, FALSE, 1, name="P0"),
	mxData(observed=demoOneFactor, type="raw"),
	imxExpectationStateSpace("A", "B", "C", "D", "Q", "R", "x0", "P0"),
	mxFitFunctionML()
)
ssRun <- mxRun(ssModel)
summary(ssRun)
# Note the Autoregressive parameter (A matrix) is near zero
#  as it should be for the independent rows of data here.

}
