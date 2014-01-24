%
%   Copyright 2007-2012 The OpenMx Project
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

\name{mxErrorPool}
\alias{mxErrorPool}

\title{Query the Error Pool}

\description{
   Retrieve models from the pool that did not complete successfully.
}

\usage{
mxErrorPool(modelnames = NA, reset = FALSE)
}

\arguments{
    \item{modelnames}{Either NA or a character vector of model names.}    
	\item{reset}{Either TRUE or FALSE.}
}

\details{
	If \sQuote{modelnames} is NA, then the list of all error models will be returned. Otherwise a subset of models will be returned, basedon the model names passed in as a argument.  If \sQuote{reset} is TRUE, then the error pool is reset to the empty list.
}

\references{
The OpenMx User's guide can be found at http://openmx.psyc.virginia.edu/documentation.
}
