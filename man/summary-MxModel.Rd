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

\name{summary-MxModel}
\alias{summary,MxModel-method}
\alias{summary}

\title{Model Summary}

\description{
   This function returns summary statistics of a model. It is usually invoked after a model has been run through the optimizer.
}


\usage{
summary(object, ...)
}

\arguments{
   \item{object}{A MxModel object.}
   \item{...}{Any number of named arguments (see below).}
}

\details{

The following named arguments are supported by the summary method:
\describe{
\item{numObs}{Numeric. Specify the total number of observations for the model.}
\item{numStats}{Numeric. Specify the total number of observed statistics for the model.}
\item{SaturatedLikelihood}{Numeric or MxModel object. Specify a saturated likelihood for testing.}
\item{indep}{Logical. Set to FALSE to ignore independent submodels in summary.}
}

}

\references{
The OpenMx User's guide can be found at http://openmx.psyc.virginia.edu/documentation.
}

\examples{

  model <- mxModel()
  modelOut <- mxRun(model)
  # compute a summary and store in variable "statistics"
  statistics <- summary(modelOut)

  # specify a saturated likelihood for testing
  summary(modelOut, SaturatedLikelihood=300)

}
