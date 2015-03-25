\name{mxExpectationGREML}
\alias{mxExpectationGREML}

\title{
Create MxExpectationGREML Object
}
\description{
This function creates a new \code{\link{MxExpectationGREML}} object.
}
\usage{
mxExpectationGREML(V, X="X", y="y", dV=character(0))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{V}{Character string; the name of the \code{\link{MxAlgebra}} or \code{\link{MxMatrix}} to serve as the 'V' matrix (the model-expected covariance matrix).}
  \item{X}{Character string; the name of the \code{\link{MxAlgebra}} or \code{\link{MxMatrix}} to serve as the 'X' matrix (the matrix of fixed covariates).}
  \item{y}{Character string; the name of the \code{\link{MxAlgebra}} or \code{\link{MxMatrix}} to serve as the 'y' matrix (the vector of phenotypic observations).}
  \item{dV}{Vector of character strings; defaults to a character vector of length zero.  If a value of non-zero length is provided, it must be a \emph{named} character vector.  This vector's names must be the labels of free parameters in the model.  The vector's elements (i.e., the character strings themselves) must be the names of \code{\link{MxAlgebra}} or \code{\link{MxMatrix}} object(s), each of which equals the first partial derivatives of the 'V' matrix with respect to the corresponding free parameter.}
}
\details{
"GREML" stands for "genomic-relatedness-matrix restricted maximum-likelihood."  In the strictest sense of the term, it refers to genetic variance-component estimation from matrices of subjects' pairwise degree of genetic relatedness, as calculated from genome-wide marker data.  It is from this original motivation that some of the terminology originates, such as calling 'y' the "phenotype" vector.  However, OpenMx's implementation of GREML is applicable for analyses from any subject-matter domain, and in which the following assumptions are reasonable:
\enumerate{
  \item Conditional on 'X' (the covariates), the phenotype vector (response variable) 'y' is a single realization from a multivariate-normal distribution having (in general) a dense covariance matrix, 'V.'
  \item The parameters of the covariance matrix, such as variance components, are of primary interest.
  \item The random effects are normally distributed.
  \item Weighted least-squares regression, using the inverse of 'V' as a weight matrix, is an adequate model for the phenotypic means.  Note that the regression coefficients are not actually free parameters to be numerically optimized.}
  
Computationally, the chief distinguishing feature of an OpenMx GREML analysis is that the phenotype vector, 'y,' is a single realization of a random vector that, in general, cannot be partitioned into independent subvectors.  For this reason, definition variables are not compatible (and should be unnecessary with) GREML expectation.  GREML expectation can still be used if the covariance matrix is sparse, but as of this writing, OpenMx does not take advantage of the sparseness to improved performance.  Because of the limitations of restricted maximum likelihood, GREML expectation is presently incompatible with ordinal variables.  Finally, as of this writing, GREML-specified models do not use \code{\link{MxData}} objects.  However, \code{\link{mxGREMLStarter}()} is a function that can produce the 'X' and 'y' matrices from a raw dataset.

It is necessary that the 'y,' 'X,' and 'V' matrices have the same number or rows, or at least, that they will have the same number of rows when incomplete observations are trimmed from 'V' at runtime; see \code{\link{mxFitFunctionGREML}()} and \code{\link{mxGREMLStarter}()} for details.  The 'X' and 'y' matrices should not contain any \code{NA}s at runtime.  The 'X' matrix should never depend upon free parameters, and generally, neither should the 'y' matrix.

Making effective use of argument \code{dV} will usually require a custom \code{\link{mxComputeSequence}()}.  The derivatives of the REML loglikelihood function with respect to parameters can be internally computed from the derivatives of the 'V' matrix supplied via \code{dV}.  These loglikelihood derivatives will be valid as long as (1) the derivatives of 'V' evaluate to symmetric matrices the same size as 'V,' (2) neither 'X' nor 'y' depend upon free parameters, and (3) the model contains no \code{\link{MxConstraint}}s.
}
\value{
Returns a new object of class \code{\link{MxExpectationGREML}}.
}
\references{
One of the first uses of the acronym "GREML":\cr
Benjamin DJ, Cesarini D, van der Loos MJHM, Dawes CT, Koellinger PD, et al. (2012) The genetic architecture of economic and political preferences. Proceedings of the National Academy of Sciences 109: 8026-8031. doi: 10.1073/pnas.1120666109

The OpenMx User's guide can be found at http://openmx.psyc.virginia.edu/documentation.
}

\seealso{
See \code{\link{MxExpectationGREML}} for the S4 class created by \code{mxExpectationGREML()}. More information about the OpenMx package may be found \link[=OpenMx]{here}. 
}
\examples{
dat <- cbind(rnorm(100),rep(1,100))
colnames(dat) <- c("y","x")

ge <- mxExpectationGREML(V="V",X="X",y="y",dV=c(ve="I"))
gff <- mxFitFunctionGREML()
plan <- mxComputeSequence(steps=list(
  mxComputeNewtonRaphson(freeSet=c("Ve"),fitfunction="fitfunction"),
  mxComputeOnce('fitfunction',
    c('fit','gradient','hessian','ihessian'),freeSet=c("Ve")),
  mxComputeStandardError(freeSet=c("Ve")),
  mxComputeReportDeriv(freeSet=c("Ve"))
))

nadat <- matrix(as.double(NA),1,1,dimnames = list("asdf","asdf"))

testmod <- mxModel(
  "GREMLtest",
  mxData(observed = nadat, type="raw", numObs=0),
  mxMatrix(type = "Full", nrow = 1, ncol=1, free=TRUE,
    values = 1, labels = "ve", lbound = 0.0001, name = "Ve"),
  mxMatrix("Iden",nrow=100,name="I",condenseSlots=TRUE),
  mxMatrix(type = "Full", nrow=100, ncol=1, free=TRUE,
    values=dat[,2], name = "X",condenseSlots=TRUE),
  mxMatrix(type = "Full", nrow=100, ncol=1, free=TRUE,
    values=dat[,1], name="y",condenseSlots=TRUE),
  mxAlgebra(I \%x\% Ve,name="V"),
  ge,
  gff,
  plan
)
str(testmod)
}
