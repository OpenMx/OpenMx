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

\name{mxExpectationStateSpaceContinuousTime}
\alias{mxExpectationStateSpaceContinuousTime}
\alias{mxExpectationSSCT}

\title{Create an MxExpectationStateSpace Object}

\description{
This function creates a new MxExpectationStateSpace object.
}

\usage{
mxExpectationStateSpaceContinuousTime(A, B, C, D, Q, R, x0, P0, u, t = NA,
                        dimnames = NA, thresholds = NA, threshnames = dimnames)
mxExpectationSSCT(A, B, C, D, Q, R, x0, P0, u, t = NA,
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
   \item{t}{A character string indicating the name of the 't' matrix.}
   \item{dimnames}{An optional character vector to be assigned to the row names of the 'C' matrix.}
   \item{thresholds}{Not Yet Implemented. An optional character string indicating the name of the thresholds matrix.}
   \item{threshnames}{Not Yet Implemented. An optional character vector to be assigned to the column names of the thresholds matrix.}
}

\details{
The \code{mxExpectationStateSpaceContinuousTime} and \code{mxExpectationSSCT} functions are identical.  The latter is simply an abbreviated name.  When using the former, tab completion is strongly encouraged to save tedious typing.  Both of these functions are wrappers for the \link{mxExpectationStateSpace} function, which could be used for both discrete and continuous time modeling.  However, there is a strong possibility of misunderstanding the model parameters when switching between discrete time and continuous time.  The expectation matrices have the same names, but mean importantly different things so caution is warranted.  The best practice is to use \code{mxExpectationStateSpace} for discrete time models, and \code{mxExpectationStateSpaceContinuousTime} for continuous time models.

Expectation functions define the way that model expectations are calculated. That is to say, expectation functions define how a set of model matrices get turned into expectations for the data.  When used in conjunction with the \link{mxFitFunctionML}, the mxExpectationStateSpace uses maximum likelihood prediction error decomposition (PED) to obtain estimates of free parameters in a model of the raw \link{MxData} object. Continuous time state space expectations treat the raw data as a multivariate time series of possibly unevenly spaced times with each row corresponding to a single occasion.  Continuous time state space expectations implement a hybrid Kalman filter to produce expectations.  The hybrid Kalman filter uses a Kalman-Bucy filter for the prediction step and the classical Kalman filter for the update step.  It is a hybrid between the classical Kalman filter used for the discrete (but possibly unequally spaced) measurement occastions and the continous time Kalman-Bucy filter for latent variable predictions.

Missing data handling is implemented in the same fashion as full information maximum likelihood for partially missing rows of data.  Additionally, completely missing rows of data are handled by only using the prediction step from the Kalman-Bucy filter and omitting the update step.

This model uses notation for the model matrices commonly found in engineering and control theory.

 The 'A', 'B', 'C', 'D', 'Q', 'R', 'x0', and 'P0' arguments must be the names of \link{MxMatrix} or \link{MxAlgebra}objects with the associated properties of the A, B, C, D, Q, R, x0, and P0 matrices in the state space modeling approach.  The 't' matrix must be a 1x1 matrix using definition variables that gives the times at which measurements occurred.

The state space expectation is defined by the following model equations.

\deqn{\frac{d}{dt} x(t) = A x(t) + B u_t + q(t)}{dx(t)/dt = A x(t) + B u[t] + q(t)}
\deqn{y_t = C x_t + D u_t + r_t}{y[t] = C x[t] + D u[t] + r[t]}
with \eqn{q(t)} and \eqn{r_t}{r[t]} both independently and identically distributed random Gaussian (normal) variables with mean zero and covariance matrices \eqn{Q} and \eqn{R}, respectively.  Subscripts or square brackets indicate discrete indices; parentheses indicate continuous indices.  The derivative of \eqn{x(t)} with respect to \eqn{t} is \eqn{\frac{d}{dt} x(t)}{dx(t)/dt}.

The first equation is called the state equation.  It describes how the latent states change over time with a first-order linear differential equation.

The second equation is called the output equation.  It describes how the latent states relate to the observed states at a single point in time.  The output equation shows how the observed output is produced by the latent states.  Also, the output equation in state space modeling is directly analogous to the measurement model in LISREL structural equation modeling.

Note that the covariates, \eqn{u}, have "instantaneous" effects on both the state and output equations.  If lagged effects are desired, then the user must create a lagged covariate by shifting their observed variable to the desired lag.

The state and output equations, together with some minimal assumptions and the Kalman filter, imply a new expected covariance matrix and means vector for every row of data.  The expected covariance matrix of row \eqn{t} is

\deqn{S_{t} = C ( A P_{t-1} A^{\sf T} + Q ) C^{\sf T} + R}{S[t] = C ( A P[t-1] A^T + Q ) C^T + R}

The expected means vector of row \eqn{t} is

\deqn{\hat{y}_{t} = C x_{t} + D u_{t}}{yhat[t] = C x[t] + D u[t]}

The 'dimnames' arguments takes an optional character vector.  % TODO: Decide how dimnames are handled.
%If this argument is not a single NA, then this vector be assigned to be the row names of the 'C' matrix and optionally to the '?' matrix, if the '?' matrix exists.

The 'A' argument refers to the \eqn{A} matrix in the State Space approach. This matrix gives the dynamics.  Entries in the diagonal give the strength of the influence of a variable's position on its slope.  Entries in the off-diagonal give the coupling strength from one variable to another.  The \eqn{A} matrix is sometimes called the state-transition model.

The 'B' argument refers to the \eqn{B} matrix in the State Space approach. This matrix consists of exogenous forces that influence the dynamics. Note that the covariate effect is contemporaneous: the covariate at time \eqn{t} has influence on the slope of the latent state also at time \eqn{t}.  A lagged effect can be created by lagged the observed variable.  The \eqn{B} matrix is sometimes called the control-input model.

The 'C' argument refers to the \eqn{C} matrix in the State Space approach. This matrix consists of contemporaneous regression coefficients from the latent variable in column \eqn{j} to the observed variable in row \eqn{i}.  This matrix is directly analogous to the factor loadings matrix in LISREL and Mplus models.  The \eqn{C} matrix is sometimes called the observation model.

The 'D' argument refers to the \eqn{D} matrix in the State Space approach. This matrix consists of contemporaneous regressive coefficients from the input (manifest covariate) variable \eqn{j} to the observed variable in row \eqn{i}. The \eqn{D} matrix is sometimes called the feedthrough or feedforward matrix.

The 'Q' argument refers to the \eqn{Q} matrix in the State Space approach. This matrix gives the covariance of the dynamic noise.  The dynamic noise can be thought of as unmeasured covariate inputs active at all times.  This matrix must be symmetric, diagonal, or zero.  As a special case, it is often diagonal. The \eqn{Q} matrix is the covariance of the process noise.  Just as in factor analysis and general structural equation modeling, the scale of the latent variables is usually set by fixing some factor loadings in the \eqn{C} matrix, or fixing some factor variances in the \eqn{Q} matrix.

The 'R' argument refers to the \eqn{R} matrix in the State Space approach. This matrix consists of residual covariances among the observed (manifest) variables. This matrix must be symmetric  As a special case, it is often diagonal.  The \eqn{R} matrix is the covariance of the observation noise.

The 'x0' argument refers to the \eqn{x_0}{x[0]} matrix in the State Space approach. This matrix consists of the column vector of the initial values for the latent variables. The state space expectation uses the \eqn{x_0}{x[0]} matrix as the starting point to recursively estimate the latent variables' values at each time.  These starting values can be difficult to pick, however, for sufficiently long time series they often do not greatly impact the estimation.

The 'P0' argument refers to the \eqn{P_0}{P[0]} matrix in the State Space approach. This matrix consists of the initial values of the covariances of the error in the initial latent variable estimates given in \eqn{x_0}{x[0]}. That is, the \eqn{P_0}{P[0]} matrix gives the covariance of \eqn{x_0 - xtrue_0}{x[0] - xtrue[0]} where \eqn{xtrue_0}{xtrue[0]} is the vector of true initial values.  \eqn{P_0}{P[0]} is a measure of the accuracy of the intial latent state estimates. The Kalman filter uses this initial covariance to recursively generated a new covariance for each time point based on the previous time point.  The Kalman filter updates this covariance so that it is as small as possible (minimum trace). Similar to the \eqn{x_0}{x[0]} matrix, these starting values are often difficult to choose.

The 'u' argument refers to the \eqn{u}{u} matrix in the State Space approach. This matrix consists of the inputs or manifest covariates of the state space expectation.  The \eqn{u}{u} matrix must be a column vector with the same number of rows as the \eqn{B} and \eqn{D} matrices have columns.  If no inputs are desired, \eqn{u} can be a zero matrix.  If time-varying inputs are desired, then they should be included as columns in the \link{MxData} object and referred to in the labels of the \eqn{u} matrix as definition variables.  There is an example of this below.

The 't' argument refers to the \eqn{t} matrix in the State Space approach.  This matrix should be 1x1 (1 row and 1 column) and not free.  The label for the element of this matrix should be 'data.YourTimeVariable'.  The 'data' part does not change, but 'YourTimeVariable' should be a name in your data set that gives the times at which measurement happened.  The units of time are up to you.  Your choice of time units will influence of the values of the parameters you estimate.  Also, recall that the model is given \eqn{x_0}{x[0]} and \eqn{P_0}{P[0]}.  These always happen at \eqn{t=0}.  So the first row of data happens some amount of time after zero.

The \link{MxMatrix} objects included as arguments may be of any type, but should have the properties described above. The mxExpectationStateSpace will not return an error for incorrect specification, but incorrect specification will likely lead to estimation problems or errors in the \link{mxRun} function.

mxExpectationStateSpaceContinuousTime evaluates with respect to an \link{MxData} object. The \link{MxData} object need not be referenced in the mxExpectationStateSpace function, but must be included in the \link{MxModel} object. mxExpectationStateSpace requires that the 'type' argument in the associated \link{MxData} object be equal to 'raw'.  Neighboring rows of the \link{MxData} object are treated as adjacent, equidistant time points increasing from the first to the last row.

To evaluate, place an mxExpectationStateSpaceContinuousTime object, the \link{mxData} object for which the expected covariance approximates, referenced \link{MxAlgebra} and \link{MxMatrix} objects, and optional \link{MxBounds} and \link{MxConstraint} objects in an \link{MxModel} object. This model may then be evaluated using the \link{mxRun} function. The results of the optimization can be found in the 'output' slot of the resulting model, and may be obtained using the \link{mxEval} function.
}

\value{
    Returns a new MxExpectationStateSpace object. mxExpectationStateSpace objects should be included with models with referenced \link{MxAlgebra}, \link{MxData} and \link{MxMatrix} objects.
}


\seealso{
\link{mxExpectationStateSpace}
}


\references{
K.J. \ifelse{latex}{\out{\r Astr\"om}}{\ifelse{html}{\out{&Aring;str&ouml;m}}{Astrom}} and R.M. Murray  (2010).  \emph{ Feedback Systems: An Introduction for Scientists and Engineers}. Princeton University Press.

J. Durbin and S.J. Koopman.  (2001).  \emph{Time Series Analysis by State Space Methods}. Oxford University Press.

R.E. Kalman  (1960).  A New Approach to Linear Filtering and Prediction Problems.  \emph{Basic Engineering, 82}, 35-45.

R.E. Kalman and R.S. Bucy  (1961).  New Results in Linear Filtering and Prediction Theory.  \emph{Transactions of the ASME, Series D, Journal of Basic Engineering, 83}, 95-108.

G. Petris  (2010).  An R Package for Dynamic Linear Models.  \emph{Journal of Statistical Software, 36}, 1-16.

The OpenMx User's guide can be found at http://openmx.psyc.virginia.edu/documentation.
}

\examples{
#------------------------------------------------------------------------------
# Example 1
# Undamped linear oscillator, i.e. a noisy sine wave.
#	Measurement error, but no dynamic error, single indicator.
#	This example works great.

#--------------------------------------
# Data Generation

require(OpenMx)

set.seed(405)
tlen <- 200
t <- seq(1.2, 50, length.out=tlen)

freqParam <- .5
initialCond <- matrix(c(2.5, 0))
x <- initialCond[1,1]*cos(freqParam*t)
plot(t, x, type='l')

measVar <- 1.5
y <- cbind(obs=x+rnorm(tlen, sd=sqrt(measVar)), tim=t)

plot(t, y[,1], type='l')

#--------------------------------------
# Model Specification

#Note: the bounds are here only to keep SLSQP from
# stepping too far off a cliff.  With the bounds in
# place, SLSQP finds the right solution.  Without
# the bounds, SLSQP goes crazy.


cdim <- list('obs', c('ksi', 'ksiDot'))

amat <- mxMatrix('Full', 2, 2, c(FALSE, TRUE, FALSE, TRUE), c(0, -.1, 1, -.2),
	name='A', lbound=-10)
bmat <- mxMatrix('Zero', 2, 1, name='B')
cmat <- mxMatrix('Full', 1, 2, FALSE, c(1, 0), name='C', dimnames=cdim)
dmat <- mxMatrix('Zero', 1, 1, name='D')
qmat <- mxMatrix('Zero', 2, 2, name='Q')
rmat <- mxMatrix('Diag', 1, 1, TRUE, .4, name='R', lbound=1e-6)
xmat <- mxMatrix('Full', 2, 1, TRUE, c(0, 0), name='x0', lbound=-10, ubound=10)
pmat <- mxMatrix('Diag', 2, 2, FALSE, 1, name='P0')
umat <- mxMatrix('Zero', 1, 1, name='u')
tmat <- mxMatrix('Full', 1, 1, name='time', labels='data.tim')

osc <- mxModel("LinearOscillator", 
	amat, bmat, cmat, dmat, qmat, rmat, xmat, pmat, umat, tmat,
	mxExpectationSSCT('A', 'B', 'C', 'D', 'Q', 'R', 'x0', 'P0', 'u', 'time'),
	mxFitFunctionML(),
	mxData(y, 'raw'))


oscr <- mxRun(osc)


#--------------------------------------
# Results Examination

summary(oscr)


(ssFreqParam <- mxEval(sqrt(-A[2,1]), oscr))
freqParam


(ssMeasVar <- mxEval(R, oscr))
measVar


dampingParam <- 0
(ssDampingParam <- mxEval(-A[2,2], oscr))
dampingParam

}
