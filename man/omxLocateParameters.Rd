%
%   Copyright 2007-2018 by the individuals mentioned in the source code history
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

\title{Get the location (model, matrix, row, column) and other info for a parameter}

\description{
    Returns a data.frame summarizing the \emph{free} parameters in a model,
	possibly filtered using \sQuote{labels}.

	For each located parameter, the label, model, matrix, row, col, value, 
	and lbound & ubound are given as a row in the dataframe.
	
	Duplicated labels return a row for each location in which they are found. 

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
   Note that \sQuote{NA} is a valid input to \sQuote{labels}.
}

\seealso{
\code{\link{omxGetParameters}}, \code{\link{omxSetParameters}}, \code{\link{omxAssignFirstParameters}} 
}

\examples{

A <- mxMatrix('Full', 2, 2, labels = c("A11", "A12", NA, "A11"), values= 1:4,
   free = TRUE, byrow = TRUE, name = 'A')

model <- mxModel(A, name = 'model')

# Request all free parameters in model
omxLocateParameters(model)

# Request free parameters "A11" and all NAs
omxLocateParameters(model, c("A11", NA))

# Works with submodel
B = mxMatrix(name = 'B', 'Full', 1, 2, labels = c("B11", "notme"),
  free = c(TRUE, FALSE), values= pi)
model <- mxModel(model, mxModel(B, name = 'subB'))

# nb: only returns free parameters ('notme' not shown)
omxLocateParameters(model)

}
