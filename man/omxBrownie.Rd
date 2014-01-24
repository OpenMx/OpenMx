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

\name{omxBrownie}
\alias{omxBrownie}

\title{Make Brownies in OpenMx}

\description{
  This function returns a brownie recipe.
}

\usage{
omxBrownie(quantity=1, walnuts=TRUE)
}

\arguments{
  \item{quantity}{Number of batches of brownies desired. Defaults to one.}
  \item{walnuts}{Logical. Indicates whether walnuts are to be included in the brownies. Defaults to TRUE.}
}

\details{
	Returns a brownie recipe. Alter the `quantity` variable to make more pans of brownies. Ingredients, equipment and procedure are listed, but neither ingredients nor equipment are provided.
}

\value{
  Returns a brownie recipe. 
}

\references{
The OpenMx User's guide can be found at http://openmx.psyc.virginia.edu/documentation.
}

\seealso{
More information about the OpenMx package may be found \link[=OpenMx]{here}. 
}

\examples{
# Return a brownie recipe
omxBrownie()
}
