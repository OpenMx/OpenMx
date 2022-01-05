%
%   Copyright 2007-2021 by the individuals mentioned in the source code history
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

\name{mxRename}
\alias{mxRename}

\title{Rename a model or submodel}

\description{
  This function re-names a model. By default, the top model will be renamed. To rename a specific model, set oldname (see examples).
  Importantly, all internal references to the old model name (e.g. in algebras) will be updated to reference the new name.
}

\usage{
mxRename(model, newname, oldname = NA)
}

\arguments{
  \item{model}{a MxModel object.}
  \item{newname}{the new name of the model.}
  \item{oldname}{the name of the target model to rename. If NA then rename top model.}
}


\value{
  Return a \link{mxModel} object with the target model renamed.
}

\references{
The OpenMx User's guide can be found at \url{https://openmx.ssri.psu.edu/documentation/}.
}


\examples{
library(OpenMx)

# Create a parent model with two submodels:
modelC <- mxModel(model= 'modelC', 
	mxModel(model= 'modelA'),
	mxModel(model= 'modelB')
)

# Rename modelC (the top model) to "model1"
model1 <- mxRename(modelC, 'model_1')

# Rename submodel "modelB" to "model_2"
model1 <- mxRename(model1, oldname = 'modelB', newname = 'model_2')

model1

}
