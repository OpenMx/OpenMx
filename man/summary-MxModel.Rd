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

\name{summary-MxModel}
\alias{summary,MxModel-method}
\alias{summary}
\alias{mxSummary}

\title{Model Summary}

\description{
   This function returns summary statistics of a model after it has been run
}


\usage{
summary(object, ..., verbose=FALSE)
}

\arguments{
   \item{object}{A MxModel object.}
   \item{...}{Any number of named arguments (see below).}
   \item{verbose}{logical. Changes the printing style for summary (see Details)}
}

\details{

mxSummary allows the user to set or override the following parameters of the model:

\describe{
\item{numObs}{Numeric. Specify the total number of observations for the model.}
\item{numStats}{Numeric. Specify the total number of observed statistics for the model.}
\item{SaturatedLikelihood}{Numeric or MxModel object. Specify a saturated likelihood for testing.}
\item{SaturatedDoF}{Numeric. When SaturatedLikelihood is numeric, specify the degrees of freedom of the saturated likelihood for testing.}
\item{indep}{Logical. Set to FALSE to ignore independent submodels in summary.}
}

The \code{verbose} argument changes the printing style for the \code{summary} of a model.  When \code{verbose=FALSE}, a relatively minimal amount of information is printed: the free parameters, the likelihood, and a few fit indices.  When more information is available, more is printed.  For example, when the model has a saturated likelihood, several additional fit indices are printed.  On the other hand, when \code{verbose=TRUE}, the compute plan, the data summary, and additional timing information are always printed.  Moreover, available fit indices are printed regarless of whether or not they are defined.  The undefined fit indices are printed as \code{NA}.  Running a saturated model and including it with the call to \code{summary} will define these fit indices and they will dislay meaningful values.  It should be noted that the \code{verbose} argument only changes the printing style, all of the same information is calculated and exists in the output of \code{summary}.  More information is displayed when \code{verbose=TRUE}, and less when \code{verbose=FALSE}.
   
This function can report Error codes as follows:
   
   \itemize{
   \item 1: The final iterate satisfies the optimality conditions to the accuracy requested, but the sequence of iterates has not yet converged. NPSOL was terminated because no further improvement could be made in the merit function (Mx status GREEN)
   \item 2: The linear constraints and bounds could not be satisfied. The problem has no feasible solution.
   \item 3: The nonlinear constraints and bounds could not be satisfied. The problem may have no feasible solution.
   \item 4: The major iteration limit was reached (Mx status BLUE).
   \item 6: The model does not satisfy the first-order optimality conditions to the required accuracy, and no improved point for the merit function could be found during the final linesearch (Mx status RED)
   \item 7: The function derivates returned by funcon or funobj appear to be incorrect.
   \item 9: An input parameter was invalid
   }

}

\references{
The OpenMx User's guide can be found at http://openmx.psyc.virginia.edu/documentation.
}

\examples{

library(OpenMx)
data(demoOneFactor)  # load the demoOneFactor dataframe
manifests <- names(demoOneFactor) # set the manifest to the 5 demo variables
latents <- c("G")  # define 1 latent variable
model <- mxModel(model="One Factor", type="RAM",
    manifestVars = manifests,
    latentVars = latents,
    mxPath(from = latents, to=manifests, labels = paste("b", 1:5, sep = "")),
    mxPath(from = manifests, arrows = 2, labels = paste("u", 1:5, sep = "")),
    mxPath(from = latents, arrows = 2, free = FALSE, values = 1.0),
    mxData(cov(demoOneFactor), type = "cov", numObs = 500)
)
model <- mxRun(model) # Run the model, returning the result into model

# Show summary of the fitted model
summary(model)

# Compute the summary and store in the variable "statistics"
statistics <- summary(model)

# Access components of the summary
statistics$parameters
statistics$SaturatedLikelihood

# Specify a saturated likelihood for testing
summary(model, SaturatedLikelihood = -3000)

# Add a CI and view it in the summary
model = mxRun(mxModel(model=model, mxCI("b5")), intervals = TRUE)
summary(model)

}
