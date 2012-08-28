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

\name{mxFIMLObjective}
\alias{mxFIMLObjective}

\title{Create MxFIMLObjective Object}

\description{
   This function creates a new MxFIMLObjective object.
}

\usage{
mxFIMLObjective(covariance, means, dimnames = NA, thresholds = NA, vector = FALSE, threshnames = dimnames)
}

\arguments{
   \item{covariance}{A character string indicating the name of the expected covariance algebra.}
   \item{means}{A character string indicating the name of the expected means algebra.}
   \item{dimnames}{An optional character vector to be assigned to the dimnames of the covariance and means algebras.}
   \item{thresholds}{An optional character string indicating the name of the thresholds matrix.}
   \item{vector}{A logical value indicating whether the objective function result is the likelihood vector.}
   \item{threshnames}{An optional character vector to be assigned to the column names of the thresholds matrix.}
}

\details{
Objective functions are functions for which free parameter values are chosen such that the value of the objective function is minimized. The mxFIMLObjective function uses full-information maximum likelihood to provide maximum likelihood estimates of free parameters in the algebra defined by the 'covariance' and 'means' arguments. The 'covariance' argument takes an \link{MxAlgebra} object, which defines the expected covariance of an associated \link{MxData} object. The 'means' argument takes an \link{MxAlgebra} object, which defines the expected means of an associated \link{MxData} object. The 'dimnames' arguments takes an optional character vector.  If this argument is not a single NA, then this vector is used to assign the dimnames of the means vector as well as the row and columns dimnames of the covariance matrix.

The 'vector' argument is either TRUE or FALSE, and determines whether the objective function returns a column vector of the likelihoods, or a single -2*(log likelihood) value.

thresholds: The name of the thresholds matrix. When needed (for modelling ordinal data), this matrix should be created using  \link{mxMatrix}(). The thresholds matrix must have as many columns as there are ordinal variables in the model, and number of rows equal to one fewer than the maximum number of levels found in the ordinal variables. The starting values  of this matrix must also be set to reasonable values. Fill each column with a set of ordered start thresholds, one for each level of this column's factor levels minus 1. These thresholds may be free if you wish them to be estimated, or fixed. The unused rows in each column, if any, can be set to any value including NA.

threshnames: A character vector consisting of the variables in the thresholds matrix, i.e., the names of ordinal variables in a model. This is necessary for OpenMx to map the thresholds matrix columns onto the variables in your data. If you set the \link{dimnames} of the columns in the thresholds matrix then threshnames is not needed.

Usage Notes: 
dimnames must be supplied where the matrices referenced by the covariance and means algebras are not themselves labeled. Failure to do so leads to an error noting that the covariance or means matrix associated with the FIML objective does not contain dimnames.

mxFIMLObjective evaluates with respect to an \link{MxData} object. The \link{MxData} object need not be referenced in the mxFIMLObjective function, but must be included in the \link{MxModel} object. mxFIMLObjective requires that the 'type' argument in the associated \link{MxData} object be equal to 'raw'. Missing values are permitted in the associated \link{MxData} object.

To evaluate, place MxFIMLObjective objects, the \link{mxData} object for which the expected covariance approximates, referenced \link{MxAlgebra} and \link{MxMatrix} objects, and optional \link{MxBounds} and \link{MxConstraint} objects in an \link{MxModel} object. This model may then be evaluated using the \link{mxRun} function. 

The results of the optimization can be reported using the \link{summary} function, or accessed directly in the 'output' slot of the resulting model (i.e., modelName@output). Components of the output may be referenced using the \link{Extract} functionality.
}

\value{
Returns a new MxFIMLObjective object. MxFIMLObjective objects should be included with models with referenced \link{MxAlgebra}, \link{MxData} and \link{MxMatrix} objects.
}

\references{
The OpenMx User's guide can be found at http://openmx.psyc.virginia.edu/documentation.
}

\examples{

A <- mxMatrix(values = 0.5, nrow = 2, ncol = 1, 
	free = TRUE, name = "A")

D <- mxMatrix(type = "Diag", values = c(0, 0.5), 
	free = c(FALSE, TRUE), nrow = 2, name = "D")
	
M <- mxMatrix(type = "Zero", nrow = 1, ncol = 2, name = "M")

expectedCov <- mxAlgebra(A \%*\% t(A) + D, "expectedCov")

objective <- mxFIMLObjective("expectedCov", "M")

model <- mxModel(A, D, expectedCov, objective)

}
