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
   \item{dataset.is.yX}{Logical; defaults to \code{FALSE}.  If \code{TRUE}, then the first column of the raw dataset is taken as-is to be the 'y' phenotype vector, and the remaining columns are taken as-is to be the 'X' matrix of covariates.  In this case, \code{\link{mxGREMLDataHandler}()} is never internally called at runtime, and all other arguments besides \code{V}, \code{casesToDropFromV}, and \code{REML} are ignored.}
   \item{casesToDropFromV}{Integer vector.  Its elements are the indices of the rows and columns of covariance matrix 'V' (and if applicable, of the rows of 'yhat') to be dropped at runtime, usually because they correspond to rows of 'y' or 'X' that contained missing observations.  By default, no cases are dropped from 'V.'  Ignored unless \code{dataset.is.yX} is \code{TRUE}.}
   \item{REML}{Logical; defaults to \code{TRUE}.  Should restricted maximum-likelihood estimation be used?}
   \item{yhat}{Character string; the name of the \code{\link{MxAlgebra}} or \code{\link{MxMatrix}} to serve as the model-expected phenotypic mean vector.  It is recommended that the algebra or matrix be a column-vector, and some of OpenMx's utility functions (e.g., \code{\link{mxAutoStart}()}) assume that to be true.  Non-empty (i.e., non-default) values should only be provided if \code{REML} is \code{FALSE}.}
}

\details{
"GREML" stands for "genomic-relatedness-matrix restricted maximum-likelihood."  In the strictest sense of the term, it refers to genetic variance-component estimation from matrices of subjects' pairwise degree of genetic relatedness, as calculated from genome-wide marker data.  It is from this original motivation that some of the terminology originates, such as calling 'y' the "phenotype" vector.  However, OpenMx's implementation of GREML is applicable for restricted maximum-likelihood analyses from any subject-matter domain, and in which the following assumptions are reasonable:
\enumerate{
  \item Conditional on 'X' (the covariates), the phenotype vector (response variable) 'y' is a single realization from a multivariate-normal distribution having (in general) a dense covariance matrix, 'V'.
  \item The parameters of the covariance matrix, such as variance components, are of primary interest.
  \item The random effects are normally distributed.
  \item Weighted least-squares regression, using the inverse of 'V' as a weight matrix, is an adequate model for the phenotypic means.  Note that the regression coefficients are not actually free parameters to be numerically optimized.
  \item No variables in the model are treated both as phenotypes in 'y' and as covariates in 'X'.}
  
Computationally, the chief distinguishing feature of an OpenMx GREML analysis is that the phenotype vector, 'y,' is a single realization of a random vector that, in general, cannot be partitioned into independent subvectors.  For this reason, definition variables are not compatible (and should be unnecessary with) GREML expectation.  GREML expectation can still be used if the covariance matrix is sparse, but as of this writing, OpenMx does not take advantage of the sparseness to improve performance.  Partly because of the limitations of restricted maximum likelihood, GREML expectation is incompatible with ordinal variables; with GREML expectation, ordinal phenotypes must be treated as though they were continuous.

When \code{REML} is \code{TRUE}, the phenotypic means are always modeled "implicitly" via weighted least-squares regression onto the covariates in 'X', using the inverse of 'V' as the weight matrix.  The phenotypic means can also be modeled "implicitly" in the same way when \code{REML} is \code{FALSE}.  However, the case of \code{REML=FALSE} also allows the user to model the phenotypic means \emph{explicitly}, as 'yhat', an arbitrary \code{\link{MxAlgebra}} or \code{\link{MxMatrix}} to serve as the model-expected phenotypic mean vector.  So long as the appropriate care is taken, use of 'yhat' allows one to enter variables into one's model both as phenotypes in 'y' and covariates in 'X'.

When \code{REML} is \code{FALSE} and the user provides \code{yhat}, the assumptions of a GREML analysis may be stated as follows:
\enumerate{
  \item Phenotype vector 'y' is equal to mean vector 'yhat' plus normally distributed random effects, at least one of which is a single realization from a multivariate-normal distribution having (in general) a dense covariance matrix, 'V'.
  \item The parameters of the covariance matrix and the mean vector are both explicit free parameters of interest.}
}
\value{
Returns a new object of class \code{\link{MxExpectationGREML}}.
}
\references{
Kirkpatrick RM, Pritikin JN, Hunter MD, & Neale MC.  (2021).  Combining structural-equation modeling with genomic-relatedness matrix restricted maximum likelihood in OpenMx.  Behavior Genetics 51: 331-342.  \doi{10.1007/s10519-020-10037-5}

The first software implementation of "GREML":\cr
Yang J, Lee SH, Goddard ME, Visscher PM.  (2011).  GCTA: a tool for genome-wide complex trait analysis.  American Journal of Human Genetics 88: 76-82.  \doi{10.1016/j.ajhg.2010.11.011}

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
