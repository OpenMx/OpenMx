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

\name{mxExpectationStateSpace}
\alias{mxExpectationStateSpace}
\alias{MxExpectationStateSpace-class}
\alias{show,MxExpectationStateSpace-method}
\alias{print,MxExpectationStateSpace-method}

\title{Create an MxExpectationStateSpace Object}

\description{
This function creates a new MxExpectationStateSpace object.
}

\usage{
mxExpectationStateSpace(A="A", B="B", C="C", D="D", Q="Q", R="R", x0="x0", P0="P0", u="u",
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
   \item{u}{A character string indicating the name of the 'u' matrix.}
   \item{dimnames}{An optional character vector to be assigned to the row names of the 'C' matrix.}
   \item{thresholds}{Not Yet Implemented. An optional character string indicating the name of the thresholds matrix.}
   \item{threshnames}{Not Yet Implemented. An optional character vector to be assigned to the column names of the thresholds matrix.}
}

\details{
Expectation functions define the way that model expectations are calculated. When used in conjunction with the \link{mxFitFunctionML}, the mxExpectationStateSpace uses maximum likelihood prediction error decomposition (PED) to obtain estimates of free parameters in a model of the raw \link{MxData} object. State space expectations treat the raw data as a multivariate time series of equally spaced times with each row corresponding to a single occasion.  This is not a model of the block Toeplitz lagged autocovariance matrix.  State space expectations implement a classical Kalman filter to produce expectations.

The following alternative filters are not yet implemented: square root Kalman filter (in Cholesky or singular value decomposition form), extended Kalman filter for linear approximations to nonlinear state space models, unscented Kalman filter for highly nonlinear state space models, Kalman-Bucy filter for continuous time modeling, hybrid Kalman filter for continuous latent time with discrete observations, and Rauch-Tung-Striebel smoother for updating forecast state estimates after a complete forward pass through the data has been made.

Missing data handling is implemented in the same fashion as full information maximum likelihood for partially missing rows of data.  Additionally, completely missing rows of data are handled by only using the prediction step from the Kalman filter and omitting the update step.

This model uses notation for the model matrices commonly found in engineering and control theory.  

 The 'A', 'B', 'C', 'D', 'Q', 'R', 'x0', and 'P0' arguments must be the names of \link{MxMatrix} or \link{MxAlgebra}objects with the associated properties of the A, B, C, D, Q, R, x0, and P0 matrices in the state space modeling approach.

The state space expectation is defined by the following model equations.

\deqn{x_{t+1} = A x_t + B u_t + q_t}{x[t+1] = A x[t] + B u[t] + q[t]}
\deqn{y_t = C x_t + D u_t + r_t}{y[t] = C x[t] + D u[t] + r[t]}
with \eqn{q_t}{q[t]} and \eqn{r_t}{r[t]} both independently and identically distributed random Gaussian (normal) variables with mean zero and covariance matrices \eqn{Q} and \eqn{R}, respectively.

The first equation is called the state equation.  It describes how the latent states change over time.  Also, the state equation in state space modeling is directly analogous to the structural model in LISREL structural equation modeling.

The second equation is called the output equation.  It describes how the latent states relate to the observed states at a single point in time.  The output equation shows how the observed output is produced by the latent states.  Also, the output equation in state space modeling is directly analogous to the measurement model in LISREL structural equation modeling.

The state and output equations, together with some minimal assumptions and the Kalman filter, imply a new expected covariance matrix and means vector for every row of data.  The expected covariance matrix of row \eqn{t+1} is

\deqn{S_{t+1} = C ( A P_{t} A^{\sf T} + Q ) C^{\sf T} + R}{S[t+1] = C ( A P[t] A^T + Q ) C^T + R}

The expected means vector of row \eqn{t+1} is

\deqn{\hat{y}_{t+1} = C x_{t+1} + D u_{t+1}}{yhat[t+1] = C x[t+1] + D u[t+1]}

The 'dimnames' arguments takes an optional character vector.  % TODO: Decide how dimnames are handled.
%If this argument is not a single NA, then this vector be assigned to be the row names of the 'C' matrix and optionally to the '?' matrix, if the '?' matrix exists.

The 'A' argument refers to the \eqn{A} matrix in the State Space approach. This matrix consists of time regressive coefficients from the latent variable in column \eqn{j} at time \eqn{t} to the latent variable in row \eqn{i} at time \eqn{t+1}.  Entries in the diagonal are autoregressive coefficients.  Entries in the off-diagonal are cross-lagged regressive coefficients.  If the \eqn{A} and \eqn{B} matrices are zero matrices, then the state space model reduces to a factor analysis.  The \eqn{A} matrix is sometimes called the state-transition model.

The 'B' argument refers to the \eqn{B} matrix in the State Space approach. This matrix consists of time regressive coefficients from the input (manifest covariate) variable \eqn{j} at time \eqn{t} to the latent variable in row \eqn{i} at time \eqn{t+1}.  The \eqn{B} matrix is sometimes called the control-input model. 

The 'C' argument refers to the \eqn{C} matrix in the State Space approach. This matrix consists of contemporaneous regression coefficients from the latent variable in column \eqn{j} to the observed variable in row \eqn{i}.  This matrix is directly analogous to the factor loadings matrix in LISREL and Mplus models.  The \eqn{C} matrix is sometimes called the observation model.

The 'D' argument refers to the \eqn{D} matrix in the State Space approach. This matrix consists of contemporaneous regressive coefficients from the input (manifest covariate) variable \eqn{j} to the observed variable in row \eqn{i}. The \eqn{D} matrix is sometimes called the feedthrough or feedforward matrix.

The 'Q' argument refers to the \eqn{Q} matrix in the State Space approach. This matrix consists of residual covariances among the latent variables.  This matrix must be symmetric.  As a special case, it is often diagonal. The \eqn{Q} matrix is the covariance of the process noise.  Just as in factor analysis and general structural equation modeling, the scale of the latent variables is usually set by fixing some factor loadings in the \eqn{C} matrix, or fixing some factor variances in the \eqn{Q} matrix.

The 'R' argument refers to the \eqn{R} matrix in the State Space approach. This matrix consists of residual covariances among the observed (manifest) variables. This matrix must be symmetric  As a special case, it is often diagonal.  The \eqn{R} matrix is the covariance of the observation noise.

The 'x0' argument refers to the \eqn{x_0}{x[0]} matrix in the State Space approach. This matrix consists of the column vector of the initial values for the latent variables. The state space expectation uses the \eqn{x_0}{x[0]} matrix as the starting point to recursively estimate the latent variables' values at each time.  These starting values can be difficult to pick, however, for sufficiently long time series they often do not greatly impact the estimation.

The 'P0' argument refers to the \eqn{P_0}{P[0]} matrix in the State Space approach. This matrix consists of the initial values of the covariances of the error in the initial latent variable estimates given in \eqn{x_0}{x[0]}. That is, the \eqn{P_0}{P[0]} matrix gives the covariance of \eqn{x_0 - xtrue_0}{x[0] - xtrue[0]} where \eqn{xtrue_0}{xtrue[0]} is the vector of true initial values.  \eqn{P_0}{P[0]} is a measure of the accuracy of the intial latent state estimates. The Kalman filter uses this initial covariance to recursively generated a new covariance for each time point based on the previous time point.  The Kalman filter updates this covariance so that it is as small as possible (minimum trace). Similar to the \eqn{x_0}{x[0]} matrix, these starting values are often difficult to choose.

The 'u' argument refers to the \eqn{u}{u} matrix in the State Space approach. This matrix consists of the inputs or manifest covariates of the state space expectation.  The \eqn{u}{u} matrix must be a column vector with the same number of rows as the \eqn{B} and \eqn{D} matrices have columns.  If no inputs are desired, \eqn{u} can be a zero matrix.  If time-varying inputs are desired, then they should be included as columns in the \link{MxData} object and referred to in the labels of the \eqn{u} matrix as definition variables.  There is an example of this below.


The \link{MxMatrix} objects included as arguments may be of any type, but should have the properties described above. The mxExpectationStateSpace will not return an error for incorrect specification, but incorrect specification will likely lead to estimation problems or errors in the \link{mxRun} function.

mxExpectationStateSpace evaluates with respect to an \link{MxData} object. The \link{MxData} object need not be referenced in the mxExpectationStateSpace function, but must be included in the \link{MxModel} object. mxExpectationStateSpace requires that the 'type' argument in the associated \link{MxData} object be equal to 'raw'.  Neighboring rows of the \link{MxData} object are treated as adjacent, equidistant time points increasing from the first to the last row.

To evaluate, place mxExpectationStateSpace objects, the \link{mxData} object for which the expected covariance approximates, referenced \link{MxAlgebra} and \link{MxMatrix} objects, and optional \link{MxBounds} and \link{MxConstraint} objects in an \link{MxModel} object. This model may then be evaluated using the \link{mxRun} function. The results of the optimization can be found in the 'output' slot of the resulting model, and may be obtained using the \link{mxEval} function..
}

\value{
    Returns a new MxExpectationStateSpace object. mxExpectationStateSpace objects should be included with models with referenced \link{MxAlgebra}, \link{MxData} and \link{MxMatrix} objects.
}


\references{
K.J. \ifelse{latex}{\out{\rAstr\"om}}{\ifelse{html}{\out{&Aring;str&ouml;m}}{Astrom}} and R.M. Murray  (2010).  \emph{ Feedback Systems: An Introduction for Scientists and Engineers}. Princeton University Press.

J. Durbin and S.J. Koopman.  (2001).  \emph{Time Series Analysis by State Space Methods}. Oxford University Press.

R.E. Kalman  (1960).  A New Approach to Linear Filtering and Prediction Problems.  \emph{Basic Engineering, 82}, 35-45.

G. Petris  (2010).  An R Package for Dynamic Linear Models.  \emph{Journal of Statistical Software, 36}, 1-16.

The OpenMx User's guide can be found at http://openmx.psyc.virginia.edu/documentation.
}

\examples{
   
# Create and fit a model using mxMatrix, mxExpectationStateSpace, and mxFitFunctionML
require(OpenMx)
data(demoOneFactor)
nvar <- ncol(demoOneFactor)
varnames <- colnames(demoOneFactor)
ssModel <- mxModel(name="State Space Manual Example",
	mxMatrix("Full", 1, 1, TRUE, .3, name="A"),
	mxMatrix("Zero", 1, 1, name="B"),
	mxMatrix("Full", nvar, 1, TRUE, .6, name="C", dimnames=list(varnames, "F1")),
	mxMatrix("Zero", nvar, 1, name="D"),
	mxMatrix("Diag", 1, 1, FALSE, 1, name="Q"),
	mxMatrix("Diag", nvar, nvar, TRUE, .2, name="R"),
	mxMatrix("Zero", 1, 1, name="x0"),
	mxMatrix("Diag", 1, 1, FALSE, 1, name="P0"),
	mxMatrix("Zero", 1, 1, name="u"),
	mxData(observed=demoOneFactor, type="raw"),
	mxExpectationStateSpace("A", "B", "C", "D", "Q", "R", "x0", "P0", "u"),
	mxFitFunctionML()
)
ssRun <- mxRun(ssModel)
summary(ssRun)
# Note the freely estimated Autoregressive parameter (A matrix)
#  is near zero as it should be for the independent rows of data
#  from the factor model.

# Create and fit a model with INPUTS using mxMatrix, mxExpectationStateSpace, and mxFitFunctionML
require(OpenMx)
data(demoOneFactor)
nvar <- ncol(demoOneFactor)
varnames <- colnames(demoOneFactor)
#demoOneFactorInputs <- cbind(demoOneFactor, V1=rep(1, nrow(demoOneFactor)))
demoOneFactorInputs <- cbind(demoOneFactor, V1=rnorm(nrow(demoOneFactor)))
ssModel <- mxModel(name="State Space Inputs Manual Example",
	mxMatrix("Full", 1, 1, TRUE, .3, name="A"),
	mxMatrix("Full", 1, 1, TRUE, values=1, name="B"),
	mxMatrix("Full", nvar, 1, TRUE, .6, name="C", dimnames=list(varnames, "F1")),
	mxMatrix("Zero", nvar, 1, name="D"),
	mxMatrix("Diag", 1, 1, FALSE, 1, name="Q"),
	mxMatrix("Diag", nvar, nvar, TRUE, .2, name="R"),
	mxMatrix("Zero", 1, 1, name="x0"),
	mxMatrix("Diag", 1, 1, FALSE, 1, name="P0"),
	mxMatrix("Full", 1, 1, FALSE, labels="data.V1", name="u"),
	mxData(observed=demoOneFactorInputs, type="raw"),
	mxExpectationStateSpace("A", "B", "C", "D", "Q", "R", "x0", "P0", u="u"),
	mxFitFunctionML()
)
ssRun <- mxRun(ssModel)
summary(ssRun)
# Note the freely estimated Autoregressive parameter (A matrix)
#  and the freely estimated Control-Input parameter (B matrix)
#  are both near zero as they should be for the independent rows of data
#  from the factor model that does not have inputs, covariates,
#  or exogenous variables.

}
