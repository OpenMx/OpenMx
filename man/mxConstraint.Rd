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

The mxConstraint function should not be used to constrain free parameters to equality with other free parameters, either by name or by their position in an \link{MxMatrix} or \link{MxAlgebra} object. Including mxConstraint objects in \link{MxModel} objects will disable standard errors and the calculation of the final Hessian, as making constraints on free parameters, and thus should be avoided when standard errors are of importance or when constraints can be made using labels, lbound or ubound arguments or algebras. Free parameters in the same \link{MxModel} may be constrained to equality by giving them the same name in their respective 'labels' matrices. Similarly, parameters may be fixed to an individual element in a \link{MxModel} object or the result of an \link{MxAlgebra} object through labeling. For example, assigning a label of ``name[1,1]`` fixes the value of a parameter at the value in first row and first column of whatever matrix or algebra is named ``name``. The mxConstraint function should be used to enforce inequalities that cannot be conveyed using other methods.
}

\value{
Returns an \link{MxConstraint} object.
}

\references{
The OpenMx User's guide can be found at http://openmx.psyc.virginia.edu/documentation.
}

\seealso{
\link{MxConstraint} for the S4 class created by mxConstraint. 
}

\examples{
#Create a constraint between MxMatrices 'A' and 'B'
constraint <- mxConstraint(A > B, name = 'AdominatesB')

# Constrain a matrix of free parameters 'K' to be equal to matrix 'limit'
K <- mxMatrix(type="Full", nrow=2, ncol=2, free=TRUE, name="K")
limit <- mxMatrix(type="Full", nrow=2, ncol=2, free=FALSE, name="limit", values=1:4)

model <- mxModel("con_test", K, limit, 
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

# Constrain two free parameters to equality using labels
equal <- mxMatrix("Full", 2, 1, free=TRUE, values=1, labels="eq", name="D")

# Constrain a matrix element in to be equal to the result of an algebra
start <- mxMatrix("Full", 1, 1, free=TRUE,  values=1, labels="param", name="F")
alg   <- mxAlgebra(log(start), name="logP")
# force the fixed parameter in matrix G to be the result of the algebra
end   <- mxMatrix("Full", 1, 1, free=FALSE, values=1, labels="logP[1,1]", name="G")

}
