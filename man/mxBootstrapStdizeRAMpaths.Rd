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
  \item{returnRaw}{Whether or not to return the raw bootstrapping results (Defaults to \code{FALSE}: returning a dataframe summarizing the results).}
}

\details{
\code{mxBootstrapStdizeRAMpaths} applies \code{\link{mxStandardizeRAMpaths}} to each bootstrap replication, thus creating a distribution of standardized estimates for each nonzero path coefficient.

The default \code{bq} (bootstrap quantiles) of c(.25, .75) correspond to a 50\% CI. This default is chosen as many more
bootstraps are required to accurately estimate more extreme quantiles. For a 95\% CI, use \code{bq=c(.025,.0975)}. 

\emph{nb}: \sQuote{bcbci} stands for \sQuote{bias-corrected bootstrap confidence interval} To learn more about bcbci and quantile methods, see Efron (1982)
and Efron and Tibshirani (1994).

\emph{note 1}: It is possible (though unlikely) that the number of nonzero paths (elements of the \code{A} and \code{S} RAM matrices) may vary among bootstrap replications.  This precludes a simple summary of the standardized paths' bootstrapping results. In this rare case, if \code{returnRaw=TRUE}, a raw list of bootstrapping results is returned, with a warning. Otherwise an error is thrown.

\emph{note 2}: \code{mxBootstrapStdizeRAMpaths} ignores sub-models. To standardize bootstrapped sub-models, run it on the sub-models directly.
}

\value{
If \code{returnRaw=FALSE} (default), it returns a dataframe containing, among other things, the standardized path coefficients as estimated from the real data, their bootstrap SEs, and the lower and upper limits of a bootstrap confidence interval.  If \code{returnRaw=TRUE}, typically, a matrix containing the raw bootstrap results is returned; this matrix has one column per non-zero path coefficient, and one row for each successfully converged bootstrap replication or, if the number of paths varies between bootstraps, a raw list of results is returned.
}

\references{
Efron B.  (1982).  \emph{The Jackknife, the Bootstrap, and Other Resampling Plans}.  Philadelphia:  Society for Industrial and Applied Mathematics.

Efron B, Tibshirani RJ.  (1994).  \emph{An Introduction to the Bootstrap}.  Boca Raton: Chapman & Hall/CRC.
}

\seealso{
	\code{\link{mxBootstrap}()}, \code{\link{mxStandardizeRAMpaths}()}, \code{\link{mxBootstrapEval}}, \code{\link{mxSummary}}
}

\examples{

\donttest{
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

tmp = mxBootstrapStdizeRAMpaths(m1_booted)
#          name label matrix row col Std.Value    Boot.SE     25.0%     75.0%
# 1  CFA.A[1,7]    NA      A  x1  F1 0.8049842 0.01583737 0.7899938 0.8124311
# 2  CFA.A[2,7]    NA      A  x2  F1 0.7935255 0.01373320 0.7865666 0.8045558
# 3  CFA.A[3,7]    NA      A  x3  F1 0.7772050 0.01629684 0.7698374 0.7907878
# 4  CFA.A[4,7]    NA      A  x4  F1 0.8248493 0.01315534 0.8150299 0.8351416
# 5  CFA.A[5,7]    NA      A  x5  F1 0.7995083 0.01479210 0.7869158 0.8057788
# 6  CFA.A[6,7]    NA      A  x6  F1 0.8126734 0.01527586 0.8012809 0.8218805
# 7  CFA.S[1,1]    NA      S  x1  x1 0.3520004 0.02546392 0.3399556 0.3759097
# 8  CFA.S[2,2]    NA      S  x2  x2 0.3703173 0.02171159 0.3526899 0.3813130
# 9  CFA.S[3,3]    NA      S  x3  x3 0.3959524 0.02529583 0.3746547 0.4073505
# 10 CFA.S[4,4]    NA      S  x4  x4 0.3196237 0.02163979 0.3025384 0.3357263
# 11 CFA.S[5,5]    NA      S  x5  x5 0.3607865 0.02364008 0.3507206 0.3807635
# 12 CFA.S[6,6]    NA      S  x6  x6 0.3395619 0.02476480 0.3245124 0.3579489
# 13 CFA.S[7,7]    NA      S  F1  F1 1.0000000 0.00000000 1.0000000 1.0000000
# 14 CFA.M[1,1]    NA      M   1  x1 2.9950397 0.08745209 2.9368758 3.0430917
# 15 CFA.M[1,2]    NA      M   1  x2 2.9775235 0.07719970 2.9109289 3.0197492
# 16 CFA.M[1,3]    NA      M   1  x3 3.0133665 0.08645522 2.9598062 3.0779683
# 17 CFA.M[1,4]    NA      M   1  x4 3.0505604 0.08210810 2.9952130 3.1103674
# 18 CFA.M[1,5]    NA      M   1  x5 2.9776983 0.07973619 2.9362410 3.0311999
# 19 CFA.M[1,6]    NA      M   1  x6 2.9830050 0.07632118 2.9360469 3.0416504

}
}
