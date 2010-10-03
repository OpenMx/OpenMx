%
%   Copyright 2007-2010 The OpenMx Project
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

\name{mxConstraint}
\alias{mxConstraint}

\title{Create MxConstraint Object}

\description{
This function creates a new \link{MxConstraint} object.
}

\usage{
mxConstraint(expression, name = NA)
}

\arguments{
   \item{expression}{An R expression of matrix operators and matrix functions.}
   \item{name}{An optional character string indicating the name of the object.}
}

\details{
The mxConstraint function defines relationships between two \link{MxAlgebra} or \link{MxMatrix} objects. They are used to affect the estimation of free parameters in the referenced objects. The constraint relation is written identically to how a \link{MxAlgebra} expression would be written.  The outermost operator in this relation must be either \sQuote{<}, \sQuote{==} or \sQuote{>}. To affect an estimation or optimization, an \link{MxConstraint} object must be included in an \link{MxModel} object with all referenced \link{MxAlgebra} and \link{MxMatrix} objects.

The mxConstraint function should not be used to constrain free parameters, either by name or by their position in an \link{MxMatrix} or \link{MxAlgebra} object. Free parameters in the same \link{MxModel} are constrained to equality by giving them the same name in their respective 'labels' matrices.
}

\value{
Returns an \link{MxConstraint} object. If used as an argument in an \link{MxModel} object, the objects referenced in the 'alg1' and 'alg2' arguments must also be included prior to optimization.
}

\references{
The OpenMx User's guide can be found at http://openmx.psyc.virginia.edu/documentation.
}

\seealso{
\link{MxConstraint} for the S4 class created by mxConstraint. \link{MxAlgebra} and \link{MxMatrix} for objects which may be entered as arguments in the 'alg' and 'alg2' arguments, and \link{mxAlgebra} and \link{mxMatrix} for the functions that create them. More information about the OpenMx package may be found \link[=OpenMx]{here}. 
}

\examples{
#Create a constraint between MxMatrices 'A' and 'B'
constraint <- mxConstraint(A > B, name = 'AdominatesB')

# Constrain a matrix of free parameters 'K' to be equal to matrix 'limit'
K <- mxMatrix(type="Full", nrow=2, ncol=2, free=TRUE, name="K")
limit <- mxMatrix(type="Full", nrow=2, ncol=2, free=FALSE, name="limit", values=1:4)

model<- mxModel("con_test", K, limit, 
	mxConstraint(K == limit, name = "Klimit_equality"), 
	mxAlgebra(min(K), name="minK"), 
	mxAlgebraObjective("minK")
)
\dontrun{
fit <- mxRun(model)
fit@matrices$K@values
}
#      [,1] [,2]
# [1,]    1    3
# [2,]    2    4

}
