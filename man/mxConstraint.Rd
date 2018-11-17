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

\name{mxConstraint}
\alias{mxConstraint}
\alias{mxConstraintFromString}

\title{Create MxConstraint Object}

\description{
This function creates a new \link{MxConstraint} object.
}

\usage{
mxConstraint(expression, name=NA, ..., jac=character(0))
}

\arguments{
   \item{expression}{The \link{MxAlgebra}-like expression representing the constraint function.}
   \item{name}{An optional character string indicating the name of the object.}
   \item{...}{Not used.  Helps OpenMx catch bad input to argument \code{expression}, and requires argument \code{jac}--meant for advanced users--to be specified by name.}
   \item{jac}{An optional character string naming the \link{MxAlgebra} or \link{MxMatrix} representing the Jacobian for the constraint function.}
}

\details{
The \code{mxConstraint()} function defines relationships between two \link{MxAlgebra} or \link{MxMatrix} objects. They are used to affect the estimation of free parameters in the referenced objects. The constraint relation is written identically to how a \link{MxAlgebra} expression would be written.  The outermost operator in this relation must be either \sQuote{<}, \sQuote{==} or \sQuote{>}. To affect an estimation or optimization, an \link{MxConstraint} object must be included in an \link{MxModel} object with all referenced \link{MxAlgebra} and \link{MxMatrix} objects.

Usage Note: Use of \code{mxConstraint()} should be avoided where it is possible to achieve the constraint by equating free parameters by label or position in an \link{MxMatrix} or \link{MxAlgebra} object. Including mxConstraints in an mxModel will disable standard errors and the calculation of the final Hessian, and thus should be avoided when standard errors are of importance. Constraints also add computational overhead. If one labels two parameters the same, the optimizer has one fewer parameter to optimize. However, if one uses mxConstraint to do the same thing, both parameters remain estimated and a Lagrangian multiplier is added to maintain the constraint. This constraint also has to have its gradients computed and the order of the Hessian grows as well.  So while both approaches should work, the mxConstraint() will take longer to do so.

Alternatives to mxConstraints include using labels, lbound or ubound arguments or algebras. Free parameters in the same \link{MxModel} may be constrained to equality by giving them the same name in their respective 'labels' matrices. Similarly, parameters may be fixed to an individual element in a \link{MxModel} object or the result of an \link{MxAlgebra} object through labeling. For example, assigning a label of ``name[1,1]`` fixes the value of a parameter at the value in first row and first column of the matrix or algebra ``name``. The mxConstraint function should be used to enforce inequalities that cannot be conveyed using other methods.

Note that constraints imposed on functions of \link[=mxMatrix]{definition variables} will only be operative for a single row of the dataset, and \emph{not} for all rows of the dataset.

Argument \code{jac} is used to provide the name of an \link{MxMatrix} or \link{MxAlgebra} that equals the matrix of first derivatives--the Jacobian--of the constraint function with respect to the free parameters.  Here, the "constraint function" refers to the constraint expression in canonical form: an arbitrary matrix expression on the left-hand side of the comparator, and a matrix of zeroes with the same dimensions on the right-hand side.  The rows of the Jacobian correspond to elements of the  matrix result of the right-hand side, in column-major order.  Each row of the Jacobian is the vector of first partial derivatives, with respect to the free parameters of the MxModel, of its corresponding element.  Each column of the Jacobian corresponds to a free parameter of the MxModel; each column must be named with the label of the corresponding free parameter.  All the \link[=mxComputeGradientDescent]{gradient-descent} optimizers are able to take advantage of user-supplied Jacobians.  %TODO: provide a concrete example

If you provide analytic Jacobians for inequality constraints with SLSQP, the inequality constraints must use the "less than" comparator.  Providing analytic Jacobians for inequality constraints with NPSOL is recommended only for expert users.  Apparently, when there are no equality constraints, NPSOL expects the Jacobians of inequality constraints to be defined as though the constraints used the "less than" comparator, whereas if there \emph{are} equality constraints, NPSOL expects the Jacobians of the inequality constraints to be defined as though the constraints used the "greater than" comparator.
}

\value{
Returns an \link{MxConstraint} object.
}

\references{
The OpenMx User's guide can be found at http://openmx.ssri.psu.edu/documentation.
}

\seealso{
\link{MxConstraint} for the S4 class created by mxConstraint. 
}

\examples{

library(OpenMx)

#Create a constraint between MxMatrices 'A' and 'B'
constraint <- mxConstraint(A > B, name = 'AdominatesB')

# Constrain matrix 'K' to be equal to matrix 'limit'

model <- mxModel(model="con_test", 
    mxMatrix(type="Full", nrow=2, ncol=2, free=TRUE, name="K"),
    mxMatrix(type="Full", nrow=2, ncol=2, free=FALSE, name="limit", values=1:4),
    mxConstraint(K == limit, name = "Klimit_equality"), 
    mxAlgebra(min(K), name="minK"), 
    mxFitFunctionAlgebra("minK")
)

fit <- mxRun(model)
fit$matrices$K$values

#      [,1] [,2]
# [1,]    1    3
# [2,]    2    4

# Constrain both free parameters of a matrix to equality using labels (both are set to "eq")
equal <- mxMatrix("Full", 2, 1, free=TRUE, values=1, labels="eq", name="D")

# Constrain a matrix element in to be equal to the result of an algebra
start <- mxMatrix("Full", 1, 1, free=TRUE,  values=1, labels="param", name="F")
alg   <- mxAlgebra(log(start), name="logP")

# Force the fixed parameter in matrix G to be the result of the algebra
end   <- mxMatrix("Full", 1, 1, free=FALSE, values=1, labels="logP[1,1]", name="G")

}
