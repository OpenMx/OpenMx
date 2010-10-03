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

\name{mxRAMObjective}
\alias{mxRAMObjective}

\title{Create MxRAMObjective Object}

\description{
   This function creates a new MxRAMObjective object.
}

\usage{
mxRAMObjective(A, S, F, M = NA, dimnames = NA, thresholds = NA, vector = FALSE)
}


\arguments{
   \item{A}{A character string indicating the name of the 'A' matrix.}
   \item{S}{A character string indicating the name of the 'S' matrix.}
   \item{F}{A character string indicating the name of the 'F' matrix.}
   \item{M}{An optional character string indicating the name of the 'M' matrix.}
   \item{dimnames}{An optional character vector to be assigned to the column names of the 'F' and 'M' matrices.}   
   \item{thresholds}{An optional character string indicating the name of the thresholds matrix.}
   \item{vector}{A logical value indicating whether the objective function result is the likelihood vector.}
}

\details{
Objective functions are functions for which free parameter values are chosen such that the value of the objective function is minimized. The mxRAMObjective provides maximum likelihood estimates of free parameters in a model of the covariance of a given \link{MxData} object. This model is defined by reticular action modeling (McArdle and McDonald, 1984). The 'A', 'S', and 'F' arguments must refer to \link{MxMatrix} objects with the associated properties of the A, S, and F matrices in the RAM modeling approach.

The 'dimnames' arguments takes an optional character vector.  If this argument is not a single NA, then this vector be assigned to be the column names of the 'F' matrix and optionally to the 'M' matrix, if the 'M' matrix exists.

The 'A' argument refers to the A or asymmetric matrix in the RAM approach. This matrix consists of all of the asymmetric paths (one-headed arrows) in the model. A free parameter in any row and column describes a regression of the variable represented by that row regressed on the variable represented in that column. 

The 'S' argument refers to the S or symmetric matrix in the RAM approach, and as such must be square. This matrix consists of all of the symmetric paths (two-headed arrows) in the model. A free parameter in any row and column describes a covariance between the variable represented by that row and the variable represented by that column. Variances are covariances between any variable at itself, which occur on the diagonal of the specified matrix.

The 'F' argument refers to the F or filter matrix in the RAM approach. If no latent variables are included in the model (i.e., the A and S matrices are of both of the same dimension as the data matrix), then the 'F' should refer to an identity matrix. If latent variables are included (i.e., the A and S matrices are not of the same dimension as the data matrix), then the 'F' argument should consist of a horizontal adhesion of an identity matrix and a matrix of zeros. 

The 'M' argument refers to the M or means matrix in the RAM approach.  It is a 1 x n matrix, where n is the number of manifest variables + the number of latent variables. The M matrix must be specified if either the mxData type is \dQuote{cov} or \dQuote{cor} and a means vector is provided, or if the mxData type is \dQuote{raw}.  Otherwise the M matrix is ignored.

The \link{MxMatrix} objects included as arguments may be of any type, but should have the properties described above. The mxRAMObjective will not return an error for incorrect specification, but incorrect specification will likely lead to estimation problems or errors in the \link{mxRun} function.

mxRAMObjective evaluates with respect to an \link{MxData} object. The \link{MxData} object need not be referenced in the mxRAMObjective function, but must be included in the \link{MxModel} object. mxRAMObjective requires that the 'type' argument in the associated \link{MxData} object be equal to 'cov', 'cor' or 'sscp'.

To evaluate, place MxRAMObjective objects, the \link{mxData} object for which the expected covariance approximates, referenced \link{MxAlgebra} and \link{MxMatrix} objects, and optional \link{MxBounds} and \link{MxConstraint} objects in an \link{MxModel} object. This model may then be evaluated using the \link{mxRun} function. The results of the optimization can be found in the 'output' slot of the resulting model, and may be obtained using the \link{mxEval} function..
}

\value{
	Returns a new MxRAMObjective object. MxRAMObjective objects should be included with models with referenced \link{MxAlgebra}, \link{MxData} and \link{MxMatrix} objects.
}

\references{
McArdle, J. J. and MacDonald, R. P. (1984). Some algebraic properties of the Reticular Action Model for moment structures. \emph{British Journal of Mathematical and Statistical Psychology, 37,} 234-251.

The OpenMx User's guide can be found at http://openmx.psyc.virginia.edu/documentation.
}

\examples{
   
   matrixA <- mxMatrix("Full", values=c(0,0.2,0,0), name="A", nrow=2, ncol=2)
   matrixS <- mxMatrix("Full", values=c(0.8,0,0,0.8), name="S", nrow=2, ncol=2, free=TRUE)
   matrixF <- mxMatrix("Full", values=c(1,0,0,1), name="F", nrow=2, ncol=2)

   # Create a RAM objective with default A, S, F matrix names
   objective <- mxRAMObjective("A", "S", "F")
   
   model <- mxModel(matrixA, matrixS, matrixF, objective)

}
