\name{mxBootstrapStdizeRAMpaths}
\alias{mxBootstrapStdizeRAMpaths}

\title{Bootstrap distribution of standardized RAM path coefficients}

\description{
Uses the distribution of a bootstrapped RAM model's raw parameters to create a bootstrapped estimate of its standardized path coefficients.

\emph{note}: Model must have already been run through \code{\link{mxBootstrap}}.
}

\usage{
mxBootstrapStdizeRAMpaths(model, bq= c(.25, .75), 
	method= c('bcbci','quantile'), returnRaw= FALSE)
}

\arguments{
  \item{model}{An MxModel that uses \link[=mxExpectationRAM]{RAM expectation} and has already been run through \code{\link{mxBootstrap}}.}
  \item{bq}{vector of 2 bootstrap quantiles corresponding to the lower and upper limits of the desired confidence interval.}
  \item{method}{One of 'bcbci' or 'quantile'.}
  \item{returnRaw}{Whether or not to return the raw bootstrapping results (Defaults to \code{FALSE}: returning a dataframe summarising the results).}
}

\details{
\code{mxBootstrapStdizeRAMpaths} applies \code{\link{mxStandardizeRAMpaths}} to each bootstrap replication, thus creating a distribution of standardized estimates for each nonzero path coefficient.

The default \code{bq} (bootstrap quantiles) of c(.25, .75) correspond to a 50% CI. This default is chosen as many more bootstraps are required to 
accurately estimate more extreme quantiles. For a 95% CI, use \code{bq=c(.025,.0975)}. 

\emph{note 1}: It is possible (though unlikely) that the number of nonzero paths (elements of the \code{A} and \code{S} RAM matrices) may vary among bootstrap replications.  This precludes a simple summary of the standardized paths' bootstrapping results. In this rare case, if \code{returnRaw=TRUE}, a raw list of bootstrapping results is returned, with a warning. Otherwise an error is thrown.

\emph{note 2}: \code{mxBootstrapStdizeRAMpaths} ignores sub-models. To standardize bootstrapped sub-models, run it on the sub-models directly.
}

\value{
If \code{returnRaw=FALSE} (default), it returns a dataframe containing, among other things, the standardized path coefficients as estimated from the real data, their bootstrap SEs, and the lower and upper limits of a bootstrap confidence interval.  If \code{returnRaw=TRUE}, typically, a matrix containing the raw bootstrap results is returned; this matrix has one column per non-zero path coefficient, and one row for each successfully converged bootstrap replication or, if the number of paths varies between bootstraps, a raw list of results is returned.
}

\seealso{
	\code{\link{mxBootstrap}()}, \code{\link{mxStandardizeRAMpaths}()}, \code{\link{mxBootstrapEval}}, \code{\link{mxSummary}}
}

\examples{

\dontrun{
require(OpenMx)
data(myFADataRaw)
manifests = c("x1","x2","x3","x4","x5","x6")

# Build and run 1-factor raw-data CFA
m1 = mxModel("CFA", type="RAM", manifestVars=manifests, latentVars="F1",
	# Factor loadings
	mxPath("F1", to = manifests, values=1),

	# Means and variances of F1 and manifests
	mxPath(from="F1", arrows=2, free=FALSE, values=1), # fix var  F1 @1
	mxPath("one", to= "F1", free= FALSE, values = 0),  # fix mean F1 @0

	# Freely-estimate means and residual variances of manifests
	mxPath(from = manifests, arrows=2, free=TRUE, values=1),
	mxPath("one", to= manifests, values = 1),

	mxData(myFADataRaw, type="raw")
)
m1 = mxRun(m1)
set.seed(170505) # Desirable for reproducibility

# ==========================
# = 1. Bootstrap the model =
# ==========================

m1_booted = mxBootstrap(m1)

# =================================================
# = 2. Estimate and accumulate a distribution of  =
# =    standardized values from each bootstrap.   =
# =================================================

tmp = mxBootstrapStdizeRAMpaths(m1)

# tmp
#                name  Estimate         SE
# 1       ind60_to_x1 1.0000000 0.00000000
# 2       ind60_to_x2 2.1803678 0.13901100
# 3       ind60_to_x3 1.8185115 0.15219019
# 4    ind60_to_dem60 1.4830002 0.39729395
# ...

}
}
