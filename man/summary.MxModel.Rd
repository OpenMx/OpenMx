%
%   Copyright 2007-2017 The OpenMx Project
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

\name{summary.MxModel}
\alias{summary.MxModel}
\alias{mxSummary}

\title{Model Summary}

\description{
   This function returns summary statistics of a model after it has been run
}


\usage{
\S3method{summary}{MxModel}(object, ..., verbose=FALSE)
}

\arguments{
   \item{object}{A MxModel object.}
   \item{...}{Any number of named arguments (see below).}
   \item{verbose}{Whether to include extra diagnostic information.}
}

\details{

mxSummary allows the user to set or override the following parameters of the model:

\describe{
\item{numObs}{Numeric. Specify the total number of observations for the model.}
\item{numStats}{Numeric. Specify the total number of observed statistics for the model.}
\item{refModels}{List of MxModel objects. Specify a saturated and independence likelihoods in single argument for testing.}
\item{SaturatedLikelihood}{Numeric or MxModel object. Specify a saturated likelihood for testing.}
\item{SaturatedDoF}{Numeric. Specify the degrees of freedom of the saturated likelihood for testing.}
\item{IndependenceLikelihood}{Numeric or MxModel object. Specify an independence likelihood for testing.}
\item{IndependenceDoF}{Numeric. Specify the degrees of freedom of the independence likelihood for testing.}
\item{indep}{Logical. Set to FALSE to ignore independent submodels in summary.}
\item{verbose}{logical. Changes the printing style for summary (see Details)}
\item{boot.quantile}{numeric. A vector of quantiles to be used to
  summarize bootstrap replication.}
}

The \code{verbose} argument changes the printing style for the \code{summary} of a model.  When \code{verbose=FALSE}, a relatively minimal amount of information is printed: the free parameters, the likelihood, and a few fit indices.  When more information is available, more is printed.  For example, when the model has a saturated likelihood, several additional fit indices are printed.  On the other hand, when \code{verbose=TRUE}, the compute plan, the data summary, and additional timing information are always printed.  Moreover, available fit indices are printed regarless of whether or not they are defined.  The undefined fit indices are printed as \code{NA}.  Running a saturated model and including it with the call to \code{summary} will define these fit indices and they will dislay meaningful values.  It should be noted that the \code{verbose} argument only changes the printing style, all of the same information is calculated and exists in the output of \code{summary}.  More information is displayed when \code{verbose=TRUE}, and less when \code{verbose=FALSE}.

The \code{Information Criteria} (AIC, BIC) are reported in a table.  The table shows different versions of the information criteria.  Each entry in the table is an AIC or BIC obtained using different penalties.  In particular, the entries of the table do not show the values of different penalties, but rather different versions of AIC and BIC.  For example the AIC is reported with both a Parameters Penalty and a Degrees of Freedom Penalty.  AIC generally takes the form \eqn{Fit + 2*k}.  With the Parameters Penalty \eqn{k} is the number of free parameters:  \eqn{AIC.param = Fit + 2*param}.  With the Degrees of Freedom Penalty, \eqn{k} is minus one times the model degrees of freedom.  So, essentially the penalty is subtracted instead of added: \eqn{AIC.param = Fit - 2*df}.  The Degrees of Freedom penalty was used in Classic Mx.  BIC is defined similarly: \eqn{Fit + k*log(N)} where \eqn{k} is either the number of free parameters or minus one times the model degrees of freedom.  The Sample-Size Adjusted BIC is only defined for the parameters penalty: \eqn{Fit + k*log((N+2)/24)}.  For raw data models, \eqn{Fit} is the minus 2 log likelihood, \eqn{-2LL}.  For covariance data, \eqn{Fit} is the Chi-squared statistic.  The \eqn{-2LL} and saturated likelihood values reported under covariance data are not necessarily meaningful on their own, but their difference yields the Chi-squared value.

The \code{refModels}, \code{SaturatedLikelihood}, \code{SaturatedDoF}, \code{IndependenceLikelihood}, and \code{IndependenceDoF} arguments can be used to obtain further fit statistics (RMSEA, CFI, TLI, Chi-Squared).  For covariance data, saturated and independence models are fitted automatically so all fit indices are reported.  For raw data, these reference models are not estimated to save computational time.  An easy way to make reference models for most cases is provided by the \link{mxRefModels} function.  When the \code{SaturatedLikelihood} or \code{IndependenceLikelihood} arguments are used, the appropriate degrees of freedom are attempted to be calculated by OpenMx.  However, depending on the model, it may sometimes behoove the user to also explicity provide the corresponding \code{SaturatedDoF} and/or \code{IndependenceDoF}.  Again, for the vast majority of cases, the \link{mxRefModels} function handles these situations effectively and conveniently.

The summary function can report Error codes as follows:
   
   \itemize{
   \item 1: The final iterate satisfies the optimality conditions to the accuracy requested, but the sequence of iterates has not yet converged. NPSOL was terminated because no further improvement could be made in the merit function (Mx status GREEN)
   \item 2: The linear constraints and bounds could not be satisfied. The problem has no feasible solution.
   \item 3: The nonlinear constraints and bounds could not be satisfied. The problem may have no feasible solution.
   \item 4: The major iteration limit was reached (Mx status BLUE).
   \item 5: The Hessian at the solution does not appear to be convex (Mx status RED).
   \item 6: The model does not satisfy the first-order optimality conditions to the required accuracy, and no improved point for the merit function could be found during the final linesearch (Mx status RED)
   \item 7: The function derivates returned by funcon or funobj appear to be incorrect.
   \item 9: An input parameter was invalid
   }

When the information matrix is available, standard errors are reported.
If the information matrix was estimated using finite differences then an
additional diagnostic column 'A' is displayed. An exclamation point in
the 'A' column indicates that the gradient appears to be asymmetric and
the standard error may not accurately reflect the variability of that
parameter. As a precaution, it is recommended that you compare the SEs
with likelihood-based confidence intervals.

For many raw data models, OpenMx does not automatically report the absolute fit indices (Chi-Squared, CFI, TLI, and RMSEA).  They are available once you fit reference models.  See the example given in \link{mxRefModels}.

With regard to RMSEA, it is important to note that OpenMx does not currently make a multigroup adjustment that some other structural equation modeling programs make.  In particular, we do not multiply the single-group RMSEA by the square root of the number of groups as suggested by Steiger (1998).  The RMSEA we use is based on the model likelihood (and degrees of freedom) as compared to the saturated model likelihood (and degrees of freedom), and we do not feel the adjustment is appropriate in this case.

OpenMx does not recommend using some fit indices.  These are GFI, AGFI, NFI, and SRMR.  The Goodness of Fit Index (GFI) and Adjusted Goodness of Fit Index (AGFI) are not recommended because they are strongly influeced by sample size and have rather high Type I error rates (Sharma, Mukherjee, Kumar, & Dillon, 2005).  The Normed Fit Index (NFI) has no penalty for model complexity.  That is, adding more parameters to a model always improves the NFI, regardless of how useful those parameters are.  Because the Non-Normed Fit Index (NNFI), also known as the Tucker-Lewis Index (TLI), does adjust for model complexity it is used instead.  Lastly, the Standardized Root Mean Square Residual (SRMR) is not reported because it (1) only applies to covariance models, having no direct extension to missing data, (2) has no penalty for model complexity, similar to the NFI, and (3) is positively biased (Hu & Bentler, 1999).

}

\seealso{
\code{\link{mxBootstrap}}
}

\references{
The OpenMx User's guide can be found at http://openmx.ssri.psu.edu/documentation.

Hu, L., & Bentler, P. M. (1999). Cutoff criteria for fit indexes in covariance structure analysis: Conventional criteria versus new alternatives. \emph{Structural Equation Modeling, 6,} 1-55.

Savalei, V. (2012). The relationship between root mean square error of approximation and model misspecification in confirmatory factor analysis models. \emph{Educational and Psychological Measurement, 72}(6), 910-932.

Sharma, S., Mukherjee, S., Kumar, A., & Dillon, W.R. (2005). A simulation study to investigate the use of cutoff values for assessing model fit in covariance structure models. \emph{Journal of Business Research, 58,} 935-43.

Steiger, J. H. (1998). A note on multiple sample extensions of the RMSEA fit index. \emph{Structural Equation Modeling: A Multidisciplinary Journal, 5(4)}, 411-419. DOI: 10.1080/10705519809540115
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
