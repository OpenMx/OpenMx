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

\name{omxLocateParameters}
\alias{omxLocateParameters}

\title{Summarize Model Parameters}

\description{
    Return a data.frame object summarizing the free parameters in the model.
}

\usage{
omxLocateParameters(model, labels = NULL, indep = FALSE)
}

\arguments{
   \item{model}{a MxModel object}
   \item{labels}{optionally specify which free parameters to retrieve.}
   \item{indep}{fetch parameters from independent submodels.}
}

\details{
   Invoking the function with the default value for the \sQuote{labels}
   argument retrieves all the free parameters.  The \sQuote{labels}
   argument can be used to select a subset of the free parameters.
   Note that \sQuote{NA} is a valid possible input to the \sQuote{labels} argument.
}

\seealso{
\code{\link{omxGetParameters}}, \code{\link{omxSetParameters}}, \code{\link{omxAssignFirstParameters}} 
}

\examples{

A <- mxMatrix('Full', 2, 2, labels = c("A11", "A12", NA, NA), values= 1:4,
   free = TRUE, byrow = TRUE, name = 'A')

model <- mxModel(A, name = 'model')

# Request all free parameters in model
omxLocateParameters(model)

# Request free parameters "A11" and all NAs
omxLocateParameters(model, c("A11", NA))


}
