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

\name{mxCI}
\alias{mxCI}

\title{Create mxCI Object}

\description{
This function creates a new \link{MxCI} object, which are used to estimate likelihood-based confidence intervals.
}

\usage{
mxCI(reference, interval = 0.95, type=c("both", "lower", "upper"))
}

\arguments{
\item{reference}{A character vector of free parameters, mxMatrices, mxMatrix elements and mxAlgebras on which confidence intervals are to be estimated, listed by name.}
\item{interval}{A scalar numeric value indicating the confidence interval to be estimated. Must be between 0 and 1. Defaults to 0.95.}
\item{type}{A character string indicating whether the upper, lower or both confidence limits are returned. Defaults to "both".}
}

\details{
The mxCI function creates \link{MxCI} objects, which can be used as arguments in \link{MxModel} objects. When models containing \link{MxCI} objects are optimized using \link{mxRun} with the \sQuote{intervals} argument set to TRUE, likelihood-based confidence intervals are returned. The likelihood-based confidence intervals calculated by \link{MxCI} objects are symmetric with respect to the change in likelihood in either direction, and are not necessarily symmetric around the parameter estimate. Estimation of confidence intervals requires both that an \link{MxCI} object be included in the model and that the \sQuote{intervals} argument of the \link{mxRun} function is set to TRUE. When estimated, confidence intervals can be accessed in the model output at \code{@output$confidenceIntervals} or by using \link{summary} on a fitted \link{MxModel} object. 

The likelihood-based confidence intervals returned using \link{MxCI} are obtained by increasing or decreasing the value of each parameter until the -2 log likelihood of the model increases by an amount corresponding to the requested interval. The confidence limit specified by the \sQuote{interval} argument is transformed into a corresponding difference in the model -2 log likelihood based on the likelihood ratio test. Thus, a requested confidence interval for a parameter will first determine the corresponding quantile from the chi-squared distribution with one degree of freedom (a value of 3.841459 when a 95 percent confidence interval is requested). That quantile will be populated into either the \sQuote{lowerdelta} slot, the \sQuote{upperdelta} slot, or both in the output \link{MxCI} object. 

Estimation of likelihood-based confidence intervals begins after optimization has been completed, with each parameter moved in the direction(s) specified in the \sQuote{type} argument until the specified increase in -2 log likelihood is reached. All other free parameters are left free for this stage of optimization. This process repeats until all confidence intervals have been calculated. The calculation of likelihood-based confidence intervals can be computationally intensive, and may add a significant amount of time to model estimation when many confidence intervals are requested.

Multiple parameters, \link[=MxMatrix]{MxMatrices} and \link[=MxAlgebra]{MxAlgebras} may be listed in the \sQuote{reference} argument. Individual elements of \link[=MxMatrix]{MxMatrices} and \link[=MxAlgebra]{MxAlgebras} may be listed as well, using the syntax \dQuote{matrix[row,col]} (see \link{Extract} for more information). Only scalar numeric values for the \sQuote{interval} argument are supported. Users requesting different confidence ranges for different parameters must use separate \link{mxCI} statements. \link{MxModel} objects can hold multiple \link{MxCI} objects, but only one confidence interval may be requested per named-entity.

Confidence interval estimation may result in model non-convergence at the confidence limit. Separate optimizer messages may be passed for each confidence limit. This has no impact on the parameter estimates themselves, but may indicate a problem with the referenced confidence limit. Model non-convergence for a particular confidence limit may indicate parameter interdependence or the influence of a parameter boundary. 

These error messages and their meanings are listed in the help for \link{mxSummary}

The validity of a confidence limit can be checked by running a model with the appropriate parameter fixed at the confidence limit in question. If the confidence limit is valid, the -2 log likelihoods of these two models should differ by the specified chi-squared criterion (as set using the \sQuote{lowerdelta} or \sQuote{upperdelta} slots in the \link{MxCI} object (you can choose which of these to set via the type parameter of mxCI). 

}
    
\value{
Returns a new \link{MxCI} object. If used as an argument in an \link{MxModel} object, the parameters, \link[=MxMatrix]{MxMatrices} and \link[=MxAlgebra]{MxAlgebras} listed in the 'reference' argument must also be included prior to optimization.
}

\references{
The OpenMx User's guide can be found at http://openmx.psyc.virginia.edu/documentation. Additional support for mxCI() can be found on the OpenMx wiki at http://openmx.psyc.virginia.edu/wiki.
}

\seealso{
\link{MxCI} for the S4 class created by mxCI. \link{MxMatrix} and \link{mxMatrix} for free parameter specification. More information about the OpenMx package may be found \link[=OpenMx]{here}. 
}

\examples{
# generate data
covariance <- matrix(c(1.0, 0.5, 0.5, 1.0), 
    nrow=2, 
    dimnames=list(c("a", "b"), c("a", "b")))
    
data <- mxData(covariance, "cov", numObs=100)

# create an expected covariance matrix
expect <- mxMatrix("Symm", 2, 2,
    free=TRUE,
    values=c(1, .5, 1),
    labels=c("var1", "cov12", "var2"),
    name="expectedCov")

# request 95 percent confidence intervals   
ci <- mxCI("Symm")

# specify the model
model <- mxModel("Confidence Interval Example",
    data, expect, ci,
    mxMLObjective("expectedCov", dimnames=c("a", "b")))

# run the model 
results <- mxRun(model, intervals=TRUE)

# view confidence intervals
print(summary(results)$CI)

# view all results
summary(results)


}
