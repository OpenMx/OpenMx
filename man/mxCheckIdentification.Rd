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

\name{mxCheckIdentification}
\alias{mxCheckIdentification}

\title{Check that a model is locally identified}

\description{
    Use the dimension of the null space of the Jacobian to determine whether or not a model is identified local to its current parameter values.
    The output is a list of the the identification status, the Jacobian, and which parameters are not identified.
}

\usage{
mxCheckIdentification(model, details=TRUE)
}

\arguments{
   \item{model}{A MxModel object or list of MxModel objects.}
   \item{details}{A MxModel object or list of MxModel objects.}
}

\details{
The mxCheckIdentification function is used to check that a model is identified.  That is, the function will tell you if the model has a unique solution in parameter space.  The function is most useful when applied to either (a) a model that has been run and had some NA standard errors, or (b) a model that has not been run but has reasonable starting values.  In the former situation, mxCheckIdentification is used as a diagnostic after a problem was indicated.  In the latter situation, mxCheckIdentification is used as a sanity check.

The method uses the Jacobian of the model expected means and the unique elements of the expected covariance matrix with respect to the free parameters.  It is the first derivative of the mapping between the free parameters and the sufficient statistics for the Normal distribution.  The method does not depend on data, but does depend on the current values of the free parameters.  Thus, it only provides local identification, not global identification.  Because the method does not depend on data, the model still could be empirically unidentified due to missing data.

The Jacobian is evaluated numerically and generally takes a few seconds, but much less than a minute.
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
latents <- c("G1", "G2")
model2 <- mxModel(model="Two Factor", type="RAM",
      manifestVars = manifests,
      latentVars = latents,
      mxPath(from = latents[1], to=manifests[1:3]),
      mxPath(from = latents[2], to=manifests[4:5]),
      mxPath(from = manifests, arrows = 2, lbound=1e-6),
      mxPath(from = latents, arrows = 2, free = FALSE, values = 1.0),
      mxData(cov(demoOneFactor), type = "cov", numObs=500)
)
fit2 <- mxRun(model2)

id2 <- mxCheckIdentification(fit2)
id2$status

model2n <- mxModel(model2, name="Non Identified Two Factor",
      mxPath(from=latents[2], arrows=2, free=TRUE, values=1)
)

mid2 <- mxCheckIdentification(model2n)

fit2n <- mxRun(model2n)

id2n <- mxCheckIdentification(fit2n)


}
