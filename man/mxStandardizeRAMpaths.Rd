\name{mxStandardizeRAMpaths}
\alias{mxStandardizeRAMpaths}
\alias{mxStandardizeRAMPaths}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Standardize RAM models' path coefficients}
\description{
Provides a dataframe containing the standardized values of all nonzero path coefficients appearing in the \code{A} and \code{S} matrices of models that use RAM expectation (either of \code{type="RAM"} or containing an explicit \code{\link{mxExpectationRAM}()} statement).  These standardized values are what the path coefficients would be if all variables in the analysis--both manifest and latent--were standardized to zero mean and unit variance.  If the means are being modeled in addition to the covariance structure, then the dataframe will also contain values of the nonzero elements of the \code{M} matrix after they have been re-scaled to standard deviation units.  Can optionally include asymptotic standard errors for the standardized and re-scaled coefficients, computed via the delta method.  Not intended for use with models that contain definition variables.
}
\usage{
mxStandardizeRAMpaths(model,SE=FALSE,cov=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{
An \code{\link{mxModel}} object, that either uses RAM expectation or contains at least one submodel that does.
}
  \item{SE}{
Logical.  Should standard errors be included with the standardized point estimates?  Defaults to \code{FALSE}.  Certain conditions are required for use of \code{SE=TRUE}; see "Details" below.
}
	\item{cov}{
	A repeated-sampling covariance matrix for the free-parameter estimates--say, from the robust "sandwich estimator," or from bootstrapping--used to calculate SEs for the standardized path coefficients.  Defaults to \code{NULL}, in which case twice the inverse of the Hessian matrix at the ML solution is used.  See below for details concerning concerning cases when \code{model} contains independent RAM submodels.
}
}
\details{
Matrix \code{A} contains the \emph{A}symmetric paths, i.e. the single-headed arrows.  Matrix \code{S} contains the \emph{S}ymmetric paths, i.e. the double-headed arrows.  The function will work even if \code{\link{mxMatrix}} objects named "A" and "S" are absent from the model, since it identifies which matrices in the model have been assigned the roles of \code{A} and \code{S} in the \code{\link{mxExpectationRAM}} statement.  Note that, in models of \code{type="RAM"}, the necessary matrices and expectation statement are automatically assembled from the \code{\link{mxPath}} objects.  If present, the \code{M} matrix will contain the means of exogenous variables and the intercepts of endogenous variables.

If \code{model} contains any submodels with \code{independent=TRUE} that use RAM expectation, \code{\link{mxStandardizeRAMpaths}()} automatically applies itself recursively over those submodels.  However, if a non-\code{NULL} matrix has been supplied for argument \code{cov}, that matrix is only used for the "container" model, and is not passed as argument to the recursive calls of the function.  To provide a covariance matrix for calculating SEs in an independent submodel, use \code{mxStandardizeRAMpaths()} directly on that submodel.

Use of \code{SE=TRUE} requires that package \code{numDeriv} be installed. It also requires that \code{model} contain no \code{\link{mxConstraint}} statements.  Finally, if \code{cov=NULL}, it requires \code{model} to have a nonempty \code{hessian} element in its output slot.  There are three common reasons why the latter condition may not be met.  First, the model may not have been run yet, i.e. it was not output by \code{\link{mxRun}()}.  Second, \code{\link{mxOption}} \code{"Hessian"} might be set to \code{"No"}.  Third, computing the Hessian matrix might possibly have been skipped per a user-defined \code{mxCompute*} statement (if any are present in the model).  If \code{model} contains RAM-expectation submodels with \code{independent=TRUE}, these conditions are checked separately for each such submodel.

In any event, using these standard errors for hypothesis-testing or forming confidence intervals is \emph{not} generally advised.  Instead, it is considered best practice to conduct likelihood-ratio tests or compute likelihood-based confidence intervals (from \code{\link{mxCI}()}), as in examples below.

The user should note that \code{mxStandardizeRAMpaths()} only cares whether an element of \code{A} or \code{S} (or \code{M}) is nonzero, and not whether it is a fixed or free parameter. So, for instance, if the function is used on a model not yet run, any free parameters in \code{A} or \code{S} initialized at zero will not appear in the function's output.

\strong{The user is warned} to interpret the output of \code{mxStandardizeRAMpaths()} cautiously if any elements of \code{A} or \code{S} depend upon "definition variables" (you have definition variables in your model if the \code{labels} of any \code{\link[=mxPath]{MxPath}} or \code{\link[=mxMatrix]{MxMatrix}} begin with \code{"data."}).  Typically, either \code{mxStandardizeRAMpaths()}'s results will be valid only for the first row of the raw dataset (and any rows identical to it), or some of the standardized coefficients will be incorrectly reported as zero.
}
\value{
If argument \code{model} is a single-group model that uses RAM expecation, then \code{mxStandardizeRAMpaths()} returns a dataframe, with one row for each nonzero path coefficient in \code{A} and \code{S} (and \code{M}, if present), and with the following columns:
  \item{name}{
Character strings that uniquely identify each nonzero path coefficient in terms of the model name, the matrix ("A", "S", or "M"), the row number, and the column number.
}
  \item{label}{
Character labels for those path coefficients that are labeled elements of an \code{\link{mxMatrix}} object, and \code{NA} for those that are not.  Note that path coefficients having the same label (and therefore the same UNstandardized value) can have different standardized values, and therefore the same label may appear more than once in this dataframe.
}
  \item{matrix}{
Character strings of "A", "S", or "M", depending on which matrix contains the given path coefficient.
}
  \item{row}{
Character.  The rownames of the matrix containing each path coefficient; row numbers are used instead if the matrix has no rownames.
}
  \item{col}{
Character.  The colnames of the matrix containing each path coefficient; column numbers are used instead if the matrix has no colnames.
}
  \item{Raw.Value}{
Numeric values of the raw (i.e., UNstandardized) path coefficients.
}
  \item{Raw.SE}{
Numeric values of the asymptotic standard errors of the raw path coefficients if if \code{SE=TRUE}, or "not_requested" otherwise. 
}
  \item{Std.Value}{
Numeric values of the standardized path coefficients.
}
  \item{Std.SE}{
Numeric values of the asymptotic standard errors of the standardized path coefficients if \code{SE=TRUE}, or "not_requested" otherwise.
}

If \code{model} is a multi-group model containing at least one submodel with RAM expectation, then \code{mxStandardizeRAMpaths()} returns a list.  The list has a number of elements equal to the number of submodels that either have RAM expectation or contain a submodel that does.  List elements corresponding to RAM-expectation submodels contain a dataframe, as described above.  List elements corresponding to "container" submodels are themselves lists, of the kind described here.
}

\seealso{
\code{\link{mxBootstrapStdizeRAMpaths}()}
}

\examples{
library(OpenMx)
data(demoOneFactor)
manifests <- names(demoOneFactor)
latents   <- c("G")
factorModel <- mxModel(model="One Factor", type="RAM",
      manifestVars = manifests,
      latentVars   = latents,
      mxPath(from=latents, to=manifests),
      mxPath(from=manifests, arrows=2, values=0.1),
      mxPath(from=latents, arrows=2,free=FALSE, values=1.0),
      mxData(cov(demoOneFactor), type="cov",numObs=500)
)
factorFit <-mxRun(factorModel)
summary(factorFit)$parameters
mxStandardizeRAMpaths(model=factorFit,SE=FALSE)

## Likelihood ratio test of variable x1's factor loading:
factorModelNull <- omxSetParameters(factorModel,labels="One Factor.A[1,6]",
                    values=0,free=FALSE)
factorFitNull <- mxRun(factorModelNull)
mxCompare(factorFit,factorFitNull)[2,"p"] #<--p-value

## Confidence intervals for all standardized paths:
factorModel2 <- mxModel(model=factorModel,
                        mxMatrix(type="Iden",nrow=nrow(factorModel$A),name="I"),
                        mxAlgebra( vec2diag(diag2vec( solve(I-A)\%*\%S\%*\%t(solve(I-A)) )\%^\%-0.5) ,
                                   name="InvSD"),
                        mxAlgebra( InvSD \%*\% A \%*\% solve(InvSD),
                                   name="Az",dimnames=dimnames(factorModel$A)),
                        mxAlgebra( InvSD \%*\% S \%*\% InvSD, 
                                   name="Sz",dimnames=dimnames(factorModel$S)),
                        mxCI(c("Az","Sz"))
)
\donttest{
factorFit2 <- mxRun(factorModel2,intervals=TRUE)
## Contains point values and confidence limits for all paths:
summary(factorFit2)$CI
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
