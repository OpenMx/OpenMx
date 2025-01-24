\name{mxExpectationGREML}
\alias{mxExpectationGREML}

\title{
Create MxExpectationGREML Object
}
\description{
This function creates a new \code{\link{MxExpectationGREML}} object.
}
\usage{
mxExpectationGREML(V, yvars=character(0), Xvars=list(), addOnes=TRUE, blockByPheno=TRUE, 
                  staggerZeroes=TRUE, dataset.is.yX=FALSE, casesToDropFromV=integer(0),
                  REML=TRUE, yhat=character(0))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{V}{Character string; the name of the \code{\link{MxAlgebra}} or \code{\link{MxMatrix}} to serve as the 'V' matrix (the model-expected covariance matrix).  Internally, the 'V' matrix is assumed to be symmetric, and its elements above the main diagonal are ignored.}
  \item{yvars, Xvars, addOnes, blockByPheno, staggerZeroes}{Passed to \code{\link{mxGREMLDataHandler}()}.}
   \item{dataset.is.yX}{Logical; defaults to \code{FALSE}.  If \code{TRUE}, then the first column of the raw dataset is taken as-is to be the 'y' phenotype vector, and the remaining columns are taken as-is to be the 'X' matrix of covariates.  In this case, \code{\link{mxGREMLDataHandler}()} is never internally called at runtime, and all other arguments besides \code{V} and \code{casesToDropFromV} are ignored.}
   \item{casesToDropFromV}{Integer vector.  Its elements are the numbers of the rows and columns of covariance matrix 'V' to be dropped at runtime, usually because they correspond to rows of 'y' or 'X' that contained missing observations.  By default, no cases are dropped from 'V.'  Ignored unless \code{dataset.is.yX=TRUE}.}
   \item{REML}{Logical; defaults to \code{TRUE}.  Should restricted maximum-likelihood estimation be used?  If \code{FALSE}, \code{yhat} must be provided, and \code{Xvars} is ignored; if \code{dataset.is.yX} is \code{TRUE}, then all data columns past the first (the phenotype vector) are ignored.}
   \item{yhat}{Character string; the name of the \code{\link{MxAlgebra}} or \code{\link{MxMatrix}} to serve as the model-expected phenotypic mean vector.  Must be provided if \code{REML} is \code{FALSE}.}
}

\details{
"GREML" stands for "genomic-relatedness-matrix restricted maximum-likelihood."  In the strictest sense of the term, it refers to genetic variance-component estimation from matrices of subjects' pairwise degree of genetic relatedness, as calculated from genome-wide marker data.  It is from this original motivation that some of the terminology originates, such as calling 'y' the "phenotype" vector.  However, OpenMx's implementation of GREML is applicable for analyses from any subject-matter domain, and in which the following assumptions are reasonable:
\enumerate{
  \item Conditional on 'X' (the covariates), the phenotype vector (response variable) 'y' is a single realization from a multivariate-normal distribution having (in general) a dense covariance matrix, 'V.'
  \item The parameters of the covariance matrix, such as variance components, are of primary interest.
  \item The random effects are normally distributed.
  \item Weighted least-squares regression, using the inverse of 'V' as a weight matrix, is an adequate model for the phenotypic means.  Note that the regression coefficients are not actually free parameters to be numerically optimized.}
  
Computationally, the chief distinguishing feature of an OpenMx GREML analysis is that the phenotype vector, 'y,' is a single realization of a random vector that, in general, cannot be partitioned into independent subvectors.  For this reason, definition variables are not compatible (and should be unnecessary with) GREML expectation.  GREML expectation can still be used if the covariance matrix is sparse, but as of this writing, OpenMx does not take advantage of the sparseness to improve performance.  Because of the limitations of restricted maximum likelihood, GREML expectation is incompatible with ordinal variables; with GREML expectation, ordinal phenotypes must be treated as though they were continuous.
}
\value{
Returns a new object of class \code{\link{MxExpectationGREML}}.
}
\references{
Kirkpatrick RM, Pritikin JN, Hunter MD, & Neale MC.  (2021).  Combining structural-equation modeling with genomic-relatedness matrix restricted maximum likelihood in OpenMx.  Behavior Genetics 51: 331-342.  \url{https://doi.org/10.1007/s10519-020-10037-5}

The first software implementation of "GREML":\cr
Yang J, Lee SH, Goddard ME, Visscher PM.  (2011).  GCTA: a tool for genome-wide complex trait analysis.  American Journal of Human Genetics 88: 76-82.  \url{https://doi.org/10.1016/j.ajhg.2010.11.011}

One of the first uses of the acronym "GREML":\cr
Benjamin DJ, Cesarini D, van der Loos MJHM, Dawes CT, Koellinger PD, et al. (2012).  The genetic architecture of economic and political preferences. Proceedings of the National Academy of Sciences 109: 8026-8031. doi: 10.1073/pnas.1120666109

The OpenMx User's guide can be found at \url{https://openmx.ssri.psu.edu/documentation/}.
}

\seealso{
See \code{\link{MxExpectationGREML}} for the S4 class created by \code{mxExpectationGREML()}. More information about the OpenMx package may be found \link[=OpenMx]{here}. 
}
\examples{
dat <- cbind(rnorm(100),rep(1,100))
colnames(dat) <- c("y","x")

ge <- mxExpectationGREML(V="V",yvars="y",Xvars=list("X"),addOnes=FALSE)
gff <- mxFitFunctionGREML(dV=c(ve="I"))
plan <- mxComputeSequence(freeSet=c("Ve"),steps=list(
  mxComputeNewtonRaphson(fitfunction="fitfunction"),
  mxComputeOnce('fitfunction',
    c('fit','gradient','hessian','ihessian')),
  mxComputeStandardError(),
  mxComputeReportDeriv(),
  mxComputeReportExpectation()
))

testmod <- mxModel(
  "GREMLtest",
  mxData(observed = dat, type="raw"),
  mxMatrix(type = "Full", nrow = 1, ncol=1, free=TRUE,
    values = 1, labels = "ve", lbound = 0.0001, name = "Ve"),
  mxMatrix("Iden",nrow=100,name="I",condenseSlots=TRUE),
  mxAlgebra(I \%x\% Ve,name="V"),
  ge,
  gff,
  plan
)
str(testmod)
}
