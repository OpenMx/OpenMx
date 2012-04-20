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

\name{mxLISRELObjective}
\alias{mxLISRELObjective}

\title{Create MxLISRELObjective Object}

\description{
   This function creates a new MxLISRELObjective object.
}

\usage{
mxLISRELObjective(LX=NA, LY=NA, BE=NA, GA=NA, PH=NA, PS=NA, TD=NA, TE=NA, TH=NA,
   TX = NA, TY = NA, KA = NA, AL = NA,
   dimnames = NA, thresholds = NA, vector = FALSE, threshnames = dimnames)
}


\arguments{
   \item{LX}{An optional character string indicating the name of the 'LX' matrix.}
   \item{LY}{An optional character string indicating the name of the 'LY' matrix.}
   \item{BE}{An optional character string indicating the name of the 'BE' matrix.}
   \item{GA}{An optional character string indicating the name of the 'GA' matrix.}
   \item{PH}{An optional character string indicating the name of the 'PH' matrix.}
   \item{PS}{An optional character string indicating the name of the 'PS' matrix.}
   \item{TD}{An optional character string indicating the name of the 'TD' matrix.}
   \item{TE}{An optional character string indicating the name of the 'TE' matrix.}
   \item{TH}{An optional character string indicating the name of the 'TH' matrix.}
   \item{TX}{An optional character string indicating the name of the 'TX' matrix.}
   \item{TY}{An optional character string indicating the name of the 'TY' matrix.}
   \item{KA}{An optional character string indicating the name of the 'KA' matrix.}
   \item{AL}{An optional character string indicating the name of the 'AL' matrix.}
   \item{dimnames}{An optional character vector that is currently ignored}%to be assigned to the column names of the 'F' and 'M' matrices.}   
   \item{thresholds}{An optional character string indicating the name of the thresholds matrix.}
   \item{vector}{A logical value indicating whether the objective function result is the likelihood vector.}
   \item{threshnames}{An optional character vector to be assigned to the column names of the thresholds matrix.}
}

\details{
Objective functions are functions for which free parameter values are chosen such that the value of the objective function is minimized. The mxLISRELObjective provides maximum likelihood estimates of free parameters in a model of the covariance of a given \link{MxData} object. This model is defined by LInear Structual RELations (LISREL; J\"{o}reskog \& S\"{o}rbom, 1982, 1996).  Arguments 'LX' through 'AL' must refer to \link{MxMatrix} objects with the associated properties of their respective matrices in the LISREL modeling approach.

The full LISREL specification has 13 matrices and is sometimes called the extended LISREL model.  It is defined by the following equations.

\deqn{\eta = \alpha + B \eta + \Gamma \xi + \zeta}
\deqn{y = \tau_y + \Lambda_y \eta + \epsilon}
\deqn{x = \tau_x + \Lambda_x \xi + \delta}

The table below is provided as a quick reference to the numerous matrices in LISREL models.  Note that NX is the number of manifest exogenous (independent) variables, the number of Xs.  NY is the number of manifest endogenous (dependent) variables, the number of Ys.  NK is the number of latent exogenous variables, the number of Ksis or Xis.  NE is the number of latent endogenous variables, the number of etas.

\tabular{clcccl}{
Matrix \tab Word \tab Abbreviation \tab Dimensions \tab Expression\tab Description \cr
\eqn{\Lambda_x} \tab Lambda x \tab LX \tab NX x NK \tab \tab Exogenous Factor Loading Matrix \cr
\eqn{\Lambda_y} \tab Lambda y \tab LY \tab NY x NE \tab \tab Endogenous Factor Loading Matrix\cr
\eqn{B} \tab BE \tab Beta \tab NE x NE \tab \tab Regressions of Latent Endogenous Variables Predicting Endogenous Variables\cr
\eqn{\Gamma} \tab Gamma \tab GA \tab NE x NK \tab \tab Regressions of Latent Exogenous Variables Predicting Endogenous Variables \cr
\eqn{\Phi} \tab Phi \tab PH \tab NK x NK \tab cov(\eqn{\xi}) \tab Covariance Matrix of Latent Exogenous Variables \cr
\eqn{\Psi} \tab Psi \tab PS \tab NE x NE \tab cov(\eqn{\zeta}) \tab Residual Covariance Matrix of Latent Endogenous Variables \cr
\eqn{\Theta_{\delta}} \tab Theta delta \tab TD \tab NX x NX \tab cov(\eqn{\delta}) \tab Residual Covariance Matrix of Manifest Exogenous Variables \cr
\eqn{\Theta_{\epsilon}} \tab Theta epsilon \tab TE \tab NY x NY \tab cov(\eqn{\epsilon}) \tab Residual Covariance Matrix of Manifest Endogenous Variables \cr
\eqn{\Theta_{\delta \epsilon}} \tab Theta delta epsilson \tab TH \tab NX x NY \tab cov(\eqn{\delta, \epsilon}) \tab Residual Covariance Matrix of Manifest Exogenous with Endogenous Variables \cr
\eqn{\tau_x} \tab tau x \tab TX \tab NX x 1 \tab \tab Residual Means of Manifest Exogenous Variables \cr
\eqn{\tau_y} \tab tau y \tab TY \tab NY x 1 \tab \tab Residual Means of Manifest Endogenous Variables \cr
\eqn{\kappa} \tab kappa \tab KA \tab NK x 1 \tab mean(\eqn{\xi}) \tab Means of Latent Exogenous Variables \cr
\eqn{\alpha} \tab alpha \tab AL \tab NE x 1 \tab \tab Residual Means of Latent Endogenous Variables \cr
}

From the extended LISREL model, several submodels can be defined.  Subtypes of the LISREL model are defined by setting some of the arguments of the LISREL objective to NA.  The first submodel is the LISREL model without means.

\deqn{\eta = B \eta + \Gamma \xi + \zeta}
\deqn{y = \Lambda_y \eta + \epsilon}
\deqn{x = \Lambda_x \xi + \delta}

The LISREL model without means requires 9 matrices: LX, LY, BE, GA, PH, PS, TD, TE, and TH.  Hence this LISREL model has TX, TY, KA, and AL as NA.  This can be accomplished be leaving these matrices at their default values.

The TX, TY, KA, and AL matrices must be specified if either the mxData type is \dQuote{cov} or \dQuote{cor} and a means vector is provided, or if the mxData type is \dQuote{raw}.  Otherwise the TX, TY, KA, and AL matrices are ignored and the model without means is estimated.

A second submodel involves only endogenous variables.

\deqn{\eta = B \eta + \zeta}
\deqn{y = \Lambda_y \eta + \epsilon}

The endogenous-only LISREL model requires 4 matrices: LY, BE, PS, and TE.  The LX, GA, PH, TD, and TH must be NA in this case.  However, means can also be specified, allowing TY and AL if the data are provided.

Another submodel involves only exogenous variables.

\deqn{x = \Lambda_x \xi + \delta}

The exogenous-model model requires 3 matrices: LX, PH, and TD.  The LY, BE, GA, PS, TE, and TH matrices must be NA.  However, means can also be specified, allowing TX and KA if the data are provided.

The model that is run depends on the matrices that are not NA.  If all 9 matrices are not NA, then the full model is run.  If only the 4 endogenous matrices are not NA, then the endogenous-only model is run.  If only the 3 exogenous matrices are not NA, then the exogenous-only model is run.  If some endogenous and exogenous matrices are not NA, but not all of them, then appropriate errors are thrown.

%The 'dimnames' arguments takes an optional character vector.  If this argument is not a single NA, then this vector be assigned to be the column names of the 'F' matrix and optionally to the 'M' matrix, if the 'M' matrix exists.

The \link{MxMatrix} objects included as arguments may be of any type, but should have the properties described above. The mxLISRELObjective will not return an error for incorrect specification, but incorrect specification will likely lead to estimation problems or errors in the \link{mxRun} function.

Like the \link{mxRAMObjective}, the mxLISRELObjective evaluates with respect to an \link{MxData} object. The \link{MxData} object need not be referenced in the mxLISRELObjective function, but must be included in the \link{MxModel} object. mxLISRELObjective requires that the 'type' argument in the associated \link{MxData} object be equal to 'cov', 'cor', 'sscp', or 'raw'.

To evaluate, place MxLISRELObjective objects, the \link{mxData} object for which the expected covariance approximates, referenced \link{MxAlgebra} and \link{MxMatrix} objects, and optional \link{MxBounds} and \link{MxConstraint} objects in an \link{MxModel} object. This model may then be evaluated using the \link{mxRun} function. The results of the optimization can be found in the 'output' slot of the resulting model, and may be obtained using the \link{mxEval} function.
}

\value{
	Returns a new MxLISRELObjective object. MxLISRELObjective objects should be included with models with referenced \link{MxAlgebra}, \link{MxData} and \link{MxMatrix} objects.
}

\references{
J\"{o}reskog, K. G. \& S\"{o}rbom, D. (1996). LISREL 8: User's Reference Guide. Lincolnwood, IL: Scientific Software International.

J\"{o}reskog, K. G. \& S\"{o}rbom, D. (1982). Recent developments in structural equation modeling. \emph{Journal of Marketing Research, 19,} 404-416.

The OpenMx User's guide can be found at http://openmx.psyc.virginia.edu/documentation.
}

\examples{
   
   #####------------------------------#####
   ##### Factor Model
   mLX <- mxMatrix("Full", values=c(.5, .6, .8, rep(0, 6), .4, .7, .5), name="LX", nrow=6, ncol=2, free=c(TRUE,TRUE,TRUE,rep(FALSE, 6),TRUE,TRUE,TRUE))
   mTD <- mxMatrix("Diag", values=c(rep(.2, 6)), name="TD", nrow=6, ncol=6, free=TRUE)
   mPH <- mxMatrix("Symm", values=c(1, .3, 1), name="PH", nrow=2, ncol=2, free=c(FALSE, TRUE, FALSE))
   
   # Create a LISREL objective with LX, TD, and PH matrix names
   objective <- mxLISRELObjective(LX="LX", TD="TD", PH="PH")
   
   model <- mxModel(mLX, mTD, mPH, objective)

}
