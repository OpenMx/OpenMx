\name{mxStandardizeLISRELpaths}
\alias{mxStandardizeLISRELpaths}
\alias{mxStandardizeLISRELPaths}
\title{Standardize LISREL models' path coefficients}
\description{
Provides a dataframe containing the standardized values of all nonzero path coefficients appearing in the 13 LISREL matrices of models that use LISREL expectation (either of \code{type="LISREL"} or containing an explicit \code{\link{mxExpectationLISREL}()} statement). These standardized values are what the path coefficients would be if all variables in the analysis--both manifest and latent, exogenous and endogenous--were standardized to zero mean and unit variance (or standard deviation units for means/intercepts). Can optionally include asymptotic standard errors for the standardized coefficients, computed via the delta method. Not intended for use with models that contain definition variables.
}
\usage{
mxStandardizeLISRELpaths(model, SE=FALSE, cov=NULL)
}
\arguments{
  \item{model}{
An \code{\link{mxModel}} object, that either uses LISREL expectation or contains at least one submodel that does.
}
  \item{SE}{
Logical. Should standard errors be included with the standardized point estimates? Defaults to \code{FALSE}. Certain conditions are required for use of \code{SE=TRUE}; see "Details" below.
}
  \item{cov}{
A repeated-sampling covariance matrix for the free-parameter estimates--say, from the robust "sandwich estimator," or from bootstrapping--used to calculate SEs for the standardized path coefficients. Defaults to \code{NULL}, in which case \code{vcov(model)} is used. See below for details concerning cases when \code{model} contains independent LISREL submodels.
}
}
\details{
The LISREL expectation uses up to 13 matrices to specify the model:
\itemize{
  \item \code{LX}: Lambda-X, factor loadings for exogenous manifest variables.
  \item \code{LY}: Lambda-Y, factor loadings for endogenous manifest variables.
  \item \code{BE}: Beta, paths between endogenous latent variables.
  \item \code{GA}: Gamma, paths from exogenous latent variables to endogenous latent variables.
  \item \code{PH}: Phi, covariance/variance matrix of exogenous latent variables.
  \item \code{PS}: Psi, covariance/variance/residual covariance matrix of endogenous latent variables.
  \item \code{TD}: Theta-Delta, residual covariance/variance matrix of exogenous manifest variables.
  \item \code{TE}: Theta-Epsilon, residual covariance/variance matrix of endogenous manifest variables.
  \item \code{TH}: Theta-Delta-Epsilon, covariance matrix between exogenous manifest residuals and endogenous manifest residuals.
  \item \code{TX}: Tau-X, intercepts of exogenous manifest variables.
  \item \code{TY}: Tau-Y, intercepts of endogenous manifest variables.
  \item \code{KA}: Kappa, means of exogenous latent variables.
  \item \code{AL}: Alpha, intercepts of endogenous latent variables.
}

The function identifies which matrices in the model have been assigned these roles in the \code{\link{mxExpectationLISREL}} statement.

Let the covariance matrix of the latent exogenous variables be \eqn{\Sigma_\xi = \Phi}, the covariance of latent endogenous variables be \eqn{\Sigma_\eta = (I - B)^{-1} (\Gamma \Phi \Gamma^T + \Psi) (I - B)^{-T}} (if exogenous latent variables are present, otherwise \eqn{\Sigma_\eta = (I - B)^{-1} \Psi (I - B)^{-T}}), and the cross-covariance between endogenous and exogenous latent variables be \eqn{\Sigma_{\eta\xi} = (I - B)^{-1} \Gamma \Phi}. Let the covariance matrices of the manifest exogenous and endogenous variables be \eqn{\Sigma_X = \Lambda_x \Phi \Lambda_x^T + \Theta_\delta} and \eqn{\Sigma_Y = \Lambda_y \Sigma_\eta \Lambda_y^T + \Theta_\epsilon}, respectively.

Let \eqn{SD_Y = \sqrt{\text{diag}(\Sigma_Y)}}, \eqn{SD_X = \sqrt{\text{diag}(\Sigma_X)}}, \eqn{SD_\eta = \sqrt{\text{diag}(\Sigma_\eta)}}, and \eqn{SD_\xi = \sqrt{\text{diag}(\Phi)}}. Let \eqn{V_Y, V_X, V_\eta, V_\xi} be the diagonal matrices of these standard deviations. The standardized matrices are then computed as:
\itemize{
  \item \eqn{\Lambda_{y, std} = V_Y^{-1} \Lambda_y V_\eta}
  \item \eqn{\Lambda_{x, std} = V_X^{-1} \Lambda_x V_\xi}
  \item \eqn{B_{std} = V_\eta^{-1} B V_\eta}
  \item \eqn{\Gamma_{std} = V_\eta^{-1} \Gamma V_\xi}
  \item \eqn{\Phi_{std} = V_\xi^{-1} \Phi V_\xi^{-1}}
  \item \eqn{\Psi_{std} = V_\eta^{-1} \Psi V_\eta^{-1}}
  \item \eqn{\Theta_{\epsilon, std} = V_Y^{-1} \Theta_\epsilon V_Y^{-1}}
  \item \eqn{\Theta_{\delta, std} = V_X^{-1} \Theta_\delta V_X^{-1}}
  \item \eqn{\Theta_{xy, std} = V_X^{-1} \Theta_{xy} V_Y^{-1}}
  \item \eqn{TX_{std} = TX V_X^{-1}}
  \item \eqn{TY_{std} = TY V_Y^{-1}}
  \item \eqn{KA_{std} = KA V_\xi^{-1}}
  \item \eqn{AL_{std} = AL V_\eta^{-1}}
}

If \code{model} contains any submodels with \code{independent=TRUE} that use LISREL expectation, \code{\link{mxStandardizeLISRELpaths}()} automatically applies itself recursively over those submodels. However, if a non-\code{NULL} matrix has been supplied for argument \code{cov}, that matrix is only used for the "container" model, and is not passed as argument to the recursive calls of the function.

Use of \code{SE=TRUE} requires that package \code{numDeriv} be installed. Further, if \code{cov=NULL}, it requires \code{vcov(model)} to return a valid result.

Using these standard errors for hypothesis-testing or forming confidence intervals is \emph{not} generally advised. Instead, it is considered best practice to conduct likelihood-ratio tests or compute likelihood-based confidence intervals (from \code{\link{mxCI}()}).

The user should note that \code{mxStandardizeLISRELpaths()} only cares whether an element of a LISREL matrix is nonzero, and not whether it is a fixed or free parameter.

\strong{The user is warned} to interpret the output of \code{mxStandardizeLISRELpaths()} cautiously if any elements depend upon "definition variables."
}
\value{
If argument \code{model} is a single-group model that uses LISREL expectation, then \code{mxStandardizeLISRELpaths()} returns a dataframe, with one row for each nonzero path coefficient in the 13 LISREL matrices, and with the following columns:
  \item{name}{
Character strings that uniquely identify each nonzero path coefficient in terms of the model name, the matrix, the row number, and the column number.
}
  \item{label}{
Character labels for those path coefficients that are labeled elements of an \code{\link{mxMatrix}} object, and \code{NA} for those that are not.
}
  \item{matrix}{
Character strings containing the name of the LISREL matrix (e.g., "LX", "PH", etc.) containing the given path coefficient.
}
  \item{row}{
Character. The rownames of the matrix containing each path coefficient; row numbers are used instead if the matrix has no rownames.
}
  \item{col}{
Character. The colnames of the matrix containing each path coefficient; column numbers are used instead if the matrix has no colnames.
}
  \item{Raw.Value}{
Numeric values of the raw (i.e., UNstandardized) path coefficients.
}
  \item{Raw.SE}{
Numeric values of the asymptotic standard errors of the raw path coefficients if \code{SE=TRUE}, or "not_requested" otherwise.
}
  \item{Std.Value}{
Numeric values of the standardized path coefficients.
}
  \item{Std.SE}{
Numeric values of the asymptotic standard errors of the standardized path coefficients if \code{SE=TRUE}, or "not_requested" otherwise.
}

If \code{model} is a multi-group model containing at least one submodel with LISREL expectation, then \code{mxStandardizeLISRELpaths()} returns a list.
}
\seealso{
\code{\link{mxStandardizeRAMpaths}}
}
\examples{
library(OpenMx)
# Set up a LISREL model
vNames <- paste("v", as.character(1:6), sep="")
dimList <- list(vNames, vNames)
covData <- matrix(
  c(0.9223099, 0.1862938, 0.4374359, 0.8959973, 0.9928430, 0.5320662,
    0.1862938, 0.2889364, 0.3927790, 0.3321639, 0.3371594, 0.4476898,
    0.4374359, 0.3927790, 1.0069552, 0.6918755, 0.7482155, 0.9013952,
    0.8959973, 0.3321639, 0.6918755, 1.8059956, 1.6142005, 0.8040448,
    0.9928430, 0.3371594, 0.7482155, 1.6142005, 1.9223567, 0.8777786,
    0.5320662, 0.4476898, 0.9013952, 0.8040448, 0.8777786, 1.3997558
    ), nrow=6, ncol=6, byrow=TRUE, dimnames=dimList)

# Create LISREL matrices
mLX <- mxMatrix("Full", values=c(.5, .6, .8, rep(0, 6), .4, .7, .5),
          name="LX", nrow=6, ncol=2,
          free=c(TRUE,TRUE,TRUE,rep(FALSE, 6),TRUE,TRUE,TRUE),
          dimnames=list(vNames, c("x1","x2")))
mTD <- mxMatrix("Diag", values=c(rep(.2, 6)), 
          name="TD", nrow=6, ncol=6, free=TRUE,
          dimnames=dimList)
mPH <- mxMatrix("Symm", values=c(1, .3, 1), 
          name="PH", nrow=2, ncol=2, free=c(FALSE, TRUE, FALSE),
          dimnames=list(c("x1","x2"),c("x1","x2")))

expFunction <- mxExpectationLISREL(LX="LX", TD="TD", PH="PH")
tmpData <- mxData(observed=covData, type="cov", numObs=100)
fitFunction <- mxFitFunctionML()
tmpModel <- mxModel(model="exampleModel",
                    mLX, mTD, mPH, expFunction, fitFunction, tmpData)
tmpModelOut <- mxRun(tmpModel)

# Standardize LISREL paths
mxStandardizeLISRELpaths(tmpModelOut, SE=TRUE)
}
