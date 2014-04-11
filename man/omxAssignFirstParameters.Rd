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

\name{omxAssignFirstParameters}
\alias{omxAssignFirstParameters}

\title{Assign First Available Values to Model Parameters}

\description{
    Sometimes you may have a free parameter with two different starting values in your model. OpenMx will not run a model until all instances of a free parameter have the same starting value. It is often sufficient to arbitrarily select one of those starting values for optimization. 

This function accomplishes that task of assigning valid starting values to the free parameters of a model. It selects an arbitrary current value (the "first" value it finds, where "first" is not defined) for each free parameter and uses that value for all instances of that parameter in the model.
}

\usage{
omxAssignFirstParameters(model, indep = FALSE)
}

\arguments{
   \item{model}{a MxModel object.}
   \item{indep}{assign parameters to independent submodels.}
}

\seealso{
\code{\link{omxGetParameters}}, \code{\link{omxSetParameters}} 
}

\examples{

A     <- mxMatrix('Full', 3, 3, values = c(1:9), labels = c('a','b', NA), free = TRUE, name = 'A')
model <- mxModel(A, name = 'model')
model <- omxAssignFirstParameters(model)

# Note: All cells with the same label now have the same start value. Note also that NAs are untouched.

model$matrices$A

# $labels
#      [,1] [,2] [,3]
# [1,] "a"  "a"  "a" 
# [2,] "b"  "b"  "b" 
# [3,] NA   NA   NA  
# 
# $values
#      [,1] [,2] [,3]
# [1,]    1    1    1
# [2,]    2    2    2
# [3,]    3    6    9

}
