%
%   Copyright 2007-2019 by the individuals mentioned in the source code history
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

\name{mxKalmanScores}
\alias{mxKalmanScores}


\title{Estimate Kalman scores and error covariance matrices}

\description{
This function creates the Kalman predicted, Kalman updated, and Rauch-Tung-Striebel smoothed latent state and error covariance estimates for an MxModel object that has an MxExpectationStateSpace object.
}

\usage{
mxKalmanScores(model, data=NA, frontend=TRUE)
}


\arguments{
   \item{model}{An MxModel object with an MxExpectationStateSpace.}
   \item{data}{An optional data.frame or matrix.}
   \item{frontend}{When TRUE, compute score in the frontend, otherwise use the backend.}
}

\details{
This is a helper function that computes the results of the classical Kalman filter.  In particular, for every row of data there is a predicted latent score, an error covariance matrix for the predicted latent scores that provides an estimate of the predictions precision, an updated latent score, and an updated error covariance matrix for the updated latent scores.  Additionally, the Rauch-Tung-Striebel (RTS) smoothed latent scores and error covariance matrices are returned.
}

\value{
    A list with components xPredicted, PPredicted, xUpdated, PUpdated, xSmoothed, PSmoothed, m2ll, and L.  When using backend scores, this list also has components for yPredicted and SPredicted which have the same number of time points as the other components but relate to the observed variables instead of the latent variables.  The rows of xPredicted, xUpdated, and xSmoothed correspond to different time points.  The columns are the different latent variables.  The third index of PPredicted, PUpdated, and PSmoothed corresponds to different times.  This works nicely with the R default print method for arrays.  At each time there is a covariance matrix of the latent variable scores.  For all items listed below, the first element goes with the zeroth time point (See example).
\describe{
	\item{xPredicted}{matrix of Kalman predicted scores}
	\item{PPredicted}{array of Kalman predicted error covariances}
	\item{xUpdated}{matrix of Kalman updated scores}
	\item{PUpdated}{array of Kalman updated error covariances}
	\item{xSmoothed}{matrix of RTS smoothed scores}
	\item{PSmoothed}{array of RTS smoothed error covariances}
	\item{m2ll}{minus 2 log likelihood}
	\item{L}{likelihood, i.e., the multivariate normal probability density}
	\item{yPredicted}{matrix of Kalman predicted scores for the observed variables, i.e., the predicted means.  Only available for backend scores.}
	\item{SPredicted}{array of Kalman predicted error covariances for the observed variables, i.e., the predicted covariances.  Only available for backend scores.}
}
}

\seealso{
\link{mxExpectationStateSpace}
}


\references{
J. Durbin and S.J. Koopman.  (2001).  \emph{Time Series Analysis by State Space Methods}. Oxford University Press.

R.E. Kalman  (1960).  A New Approach to Linear Filtering and Prediction Problems.  \emph{Basic Engineering, 82}, 35-45.

H.E. Rauch, F. Tung, C.T. Striebel.  (1965).  Maximum Likelihood Estimates of Linear Dynamic Systems.   \emph{American Institute of Aeronautics and Astronautics Journal, 3}, 1445-1450.

The OpenMx User's guide can be found at http://openmx.ssri.psu.edu/documentation.
}

\examples{
   
# Create and fit a model using mxMatrix, mxExpectationStateSpace, and mxFitFunctionML
require(OpenMx)
data(demoOneFactor)
# Use only first 50 rows, for speed of example
data <- demoOneFactor[1:50,]
nvar <- ncol(demoOneFactor)
varnames <- colnames(demoOneFactor)
ssModel <- mxModel(model="State Space Manual Example",
    mxMatrix("Full", 1, 1, TRUE, .3, name="A"),
    mxMatrix("Zero", 1, 1, name="B"),
    mxMatrix("Full", nvar, 1, TRUE, .6, name="C", dimnames=list(varnames, "F1")),
    mxMatrix("Zero", nvar, 1, name="D"),
    mxMatrix("Diag", 1, 1, FALSE, 1, name="Q"),
    mxMatrix("Diag", nvar, nvar, TRUE, .2, name="R"),
    mxMatrix("Zero", 1, 1, name="x0"),
    mxMatrix("Diag", 1, 1, FALSE, 1, name="P0"),
    mxMatrix("Zero", 1, 1, name="u"),
    mxData(observed=data, type="raw"),
    mxExpectationStateSpace("A", "B", "C", "D", "Q", "R", "x0", "P0", "u"),
    mxFitFunctionML()
)
ssRun <- mxRun(ssModel)
summary(ssRun)
# Note the freely estimated Autoregressive parameter (A matrix)
#  is near zero as it should be for the independent rows of data
#  from the factor model.

ssScores <- mxKalmanScores(ssRun)

cor(cbind(ssScores$xPredicted[,1], ssScores$xUpdated[,1], ssScores$xSmoothed[,1]))
# Because the autoregressive dynamics are near zero, the predicted and updated scores
# correlate minimally, and the updated and smoothed latent state estimates
# are extremely close.

# The first few latent predicted scores
head(ssScores$xPredicted)

# The predicted latent score for time 10
ssScores$xPredicted[10+1,]

# The error covariance of the predicted score at time 10
ssScores$PPredicted[,,10+1]

}
