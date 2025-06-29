%
%   Copyright 2007-2019 by the individuals mentioned in the source code history
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

\name{mxCheckIdentification}
\alias{mxCheckIdentification}

\title{Check that a model is locally identified}

\description{
    Use the dimension of the null space of the Jacobian to determine whether or not a model is identified local to its current parameter values.
    The output is a list of the the identification status, the Jacobian, and which parameters are not identified.  May give inaccurate results with
    \link[=mxExpectationGREML]{GREML expectations} (see "Details" below).
}

\usage{
mxCheckIdentification(model, details=TRUE, nrows=2, exhaustive=FALSE, silent=FALSE)
}

\arguments{
   \item{model}{A MxModel object or list of MxModel objects.}
   \item{details}{logical.}
   \item{nrows}{numeric. The maximum number of unique definition variable rows used to evaluate identification.}
   \item{exhaustive}{logical. Whether to use all of nrows (TRUE), or to
   stop evaluating new rows once the model is identified (FALSE).}
   \item{silent}{logical. Whether or not to print helpful-but-long
   messages.  silent=TRUE does not print messages.}
}

\details{
The mxCheckIdentification function is used to check that a model is identified.  That is, the function will tell you if the model has a unique solution in parameter space.  The function is most useful when applied to either (a) a model that has been run and had some NA standard errors, or (b) a model that has not been run but has reasonable starting values.  In the former situation, mxCheckIdentification is used as a diagnostic after a problem was indicated.  In the latter situation, mxCheckIdentification is used as a sanity check.

The method uses the Jacobian of the model expected means and the unique elements of the expected covariance matrix with respect to the free parameters.  It is the first derivative of the mapping between the free parameters and the sufficient statistics for the Normal distribution.  The method does not depend on data, but does depend on the current values of the free parameters.  Thus, it only provides local identification, not global identification.  \emph{You might get different answers about model identification depending on the free parameter values}.  Because the method does not depend on data, the model still could be empirically unidentified due to missing data.

The Jacobian is evaluated numerically and generally takes a few seconds, but much less than a minute.

Model identification should be accurate for models with linear or
nonlinear equality and inequality constraints.  When there are
constraints, mxCheckIdentification uses the Jacobian of the summary
statistics with respect to the free parameters and the Jacobian of the
summary statistics with respect to the constraints.  The combined
extended Jacobian must have rank equal to the number of free parameters
for the model to be identified.  So, a model can be identified with
constraints that is not identified without constraints.

Model identification should also be accurate for multiple group models
and models with definition variables.  In the case of multiple groups,
the Jacobian is computed for each group and they are concatenated:
summary statistics for each group go down the rows of the Jacobian; the
single set of free parameters go across the columns of the Jacobian.  In
the case of definition variables, a new set of summary statistics is
computed for each unique set of definition variable values; thus
creating a kind of pseudo-multiple group model for each set of unique
definition variable values.

For models using definition variables, the additional arguments 'nrows'
and 'exhaustive' provide instruction for how to search for
identification.  In principle, full model identification for models with
definition variables requires evaluating the Jacobian of the mapping
between the free parameters and the summary statistics \emph{for every
	unique combination of definition variable values}.  This evaluation
can take a long time.  However, in practice, the vast majority of models
are identified after evaluating one or two definition variable rows.
The current defaults attempt to capitalize on this practical situation.
By default 'mxCheckIdentification' will evaluate up to 2 unique definition
variable rows, but will stop once the model is identified.  To keep
evaluating unique definition variable values until the model is
identified or you run out of rows, set \code{nrows=Inf}.  To evaluate
all unique definition variable values, set \code{nrows=Inf} and
\code{exhaustive=TRUE}.

When TRUE, the 'details' argument provides the names of the non-identified parameters.  Otherwise, only the status and Jacobian are returned.

For most models with a \link[=mxExpectationGREML]{GREML expectation}, the function will run slowly, and use of \code{detail=TRUE} is therefore not recommended.  Additionally, for GREML \link[=mxModel]{MxModels} that use an "implicit" (i.e., "specified in terms of covariates via \code{\link{mxExpectationGREML}()} argument \code{Xvars}") model for the phenotypic mean, there are forms of unidentification that the function cannot (yet) detect, and it therefore may give false-positive reports of local identification.
}

\value{
A named list with components
\describe{
\item{status}{logical. TRUE if the model is locally identified; otherwise FALSE.}
\item{jacobian}{matrix. The numerically evaluated Jacobian.}
\item{non_identified_parameters}{vector. The free parameter names that are not identified}
}
}

\seealso{
\code{\link{mxModel}}
}

\references{
Bekker, P.A., Merckens, A., Wansbeek, T.J. (1994).  Identification, Equivalent Models and Computer Algebra.  Academic Press: Orlando, FL.

Bollen, K. A. & Bauldry, S. (2010).  Model Identification and Computer Algebra.  Sociological Methods & Research, 39, p. 127-156.
}

\examples{

require(OpenMx)

data(demoOneFactor)
manifests <- names(demoOneFactor)
latents <- "G1"
model2 <- mxModel(model="One Factor", type="RAM",
      manifestVars = manifests,
      latentVars = latents,
      mxPath(from = latents[1], to=manifests[1:5]),
      mxPath(from = manifests, arrows = 2, lbound=1e-6),
      mxPath(from = latents, arrows = 2, free = FALSE, values = 1.0),
      mxData(cov(demoOneFactor), type = "cov", numObs=500)
)
fit2 <- mxRun(model2)

id2 <- mxCheckIdentification(fit2)
id2$status
# The model is locally identified

# Build a model from the solution of the previous one
#  but now the factor variance is also free
model2n <- mxModel(fit2, name="Non Identified Two Factor",
      mxPath(from=latents[1], arrows=2, free=TRUE, values=1)
)

mid2 <- mxCheckIdentification(model2n)
mid2$non_identified_parameters
# The factor loadings and factor variance
#  are not identified.

}
