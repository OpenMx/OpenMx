#
#   Copyright 2012 The OpenMx Project
#
#   Licensed under the Apache License, Version 2.0 (the "License");
#   you may not use this file except in compliance with the License.
#   You may obtain a copy of the License at
# 
#        http://www.apache.org/licenses/LICENSE-2.0
# 
#   Unless required by applicable law or agreed to in writing, software
#   distributed under the License is distributed on an "AS IS" BASIS,
#   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#   See the License for the specific language governing permissions and
#   limitations under the License.


setClass(Class = "MxExpectationBA81",
         representation = representation(
	   ItemSpec = "list",
	   item = "MxCharOrNumber",
	   EstepItem = "MxOptionalMatrix",
	   qpoints = "numeric",
	   qwidth = "numeric",
	   mean = "MxOptionalCharOrNumber",
	   cov = "MxOptionalCharOrNumber",
	     debugInternal="logical",
	   dims = "character",
	   verbose = "integer",
	     minItemsPerScore = "integer",
	     weightColumn = "MxCharOrNumber",
	 .detectIndependence = "logical"),
         contains = "MxBaseExpectation")

setMethod("initialize", "MxExpectationBA81",
          function(.Object, ItemSpec, item, EstepItem,
		   qpoints, qwidth, mean, cov, verbose, debugInternal,
		   minItemsPerScore, weightColumn, name = 'expectation') {
            .Object@name <- name
	    .Object@ItemSpec <- ItemSpec
	    .Object@item <- item
	    .Object@EstepItem <- EstepItem
            .Object@qpoints <- qpoints
            .Object@qwidth <- qwidth
            .Object@data <- as.integer(NA)
	    .Object@mean <- mean
	    .Object@cov <- cov
	    .Object@verbose <- verbose
	    .Object@debugInternal <- debugInternal
	    .Object@minItemsPerScore <- minItemsPerScore
	    .Object@weightColumn <- weightColumn
	    .Object@.detectIndependence <- TRUE
            return(.Object)
          }
)

setMethod("genericExpDependencies", signature("MxExpectationBA81"),
	  function(.Object, dependencies) {
		  sources <- c(.Object@mean, .Object@cov, .Object@item)
		  dependencies <- imxAddDependency(sources, .Object@name, dependencies)
		  return(dependencies)
	  })

setMethod("genericExpFunConvert", signature("MxExpectationBA81"), 
	  function(.Object, flatModel, model, labelsData, dependencies) {
		  modelname <- imxReverseIdentifier(model, .Object@name)[[1]]
		  if(is.na(.Object@data)) {
			  msg <- paste(typeof(.Object),
				       "does not have a dataset associated with it in model",
				       omxQuotes(modelname))
			  stop(msg, call.=FALSE)
		  }

		  verifyMvnNames(.Object@cov, .Object@mean, "prior", flatModel, model@name, class(.Object))

		  fnames <- dimnames(flatModel[[.Object@cov]])[[1]]

		  item <- flatModel@matrices[[.Object@item]]
		  if (length(fnames) && (is.null(rownames(item)) || any(rownames(item)[1:length(fnames)] != fnames))) {
			  msg <- paste("The first", length(fnames), "rownames of item",
				       "must be", omxQuotes(fnames))
			  stop(msg)
		  }

		  name <- .Object@name
		  for (s in c("data", "item")) {
			  if (is.null(slot(.Object, s))) next
			  slot(.Object, s) <-
			    imxLocateIndex(flatModel, slot(.Object, s), name)
		  }
		  for (s in c("mean", "cov")) {
			  slot(.Object, s) <- LocateOptionalMatrix(flatModel, slot(.Object, s), name)
		  }

		  mxData <- flatModel@datasets[[.Object@data + 1]]
		  .Object@dims <- colnames(mxData@observed) # for summary

		  if (is.null(colnames(item))) {
			  stop(paste(class(.Object),
				     ": the item matrix must have column names",
				     "to match against the observed data column names."))
		  }
		  
		  dc <- match(colnames(item), .Object@dims) - 1L
		  if (any(is.na(dc))) {
			  msg <- paste("Some items are not found among the observed data columns:",
				       omxQuotes(colnames(item)[is.na(dc)]))
			  stop(msg, call.=FALSE)
		  }

		  .Object@dataColumns <- dc
		  if (!is.na(.Object@weightColumn)) {
			  wc <- match(.Object@weightColumn, .Object@dims) - 1L
			  if (is.na(wc)) {
				  msg <- paste("Weight column",.Object@weightColumn,
					       "not found in observed data columns")
				  stop(msg, call.=FALSE)
			  }
			  numObs <- sum(mxData@observed[, .Object@weightColumn])
			  if (mxData@numObs != numObs) {
				  stop(paste("For mxData, the number of observations must be set to ", numObs,
					     "; try mxData(observed=obs, type='raw',",
					     " numObs=sum(obs$", .Object@weightColumn, "), ...)", sep=""))
			  }
			  .Object@weightColumn <- wc
		  }

		  return(.Object)
	  })

setMethod("qualifyNames", signature("MxExpectationBA81"), 
	function(.Object, modelname, namespace) {
		.Object@name <- imxIdentifier(modelname, .Object@name)
		for (s in c("item", "mean", "cov")) {
			if (is.null(slot(.Object, s))) next;
			slot(.Object, s) <-
			  imxConvertIdentifier(slot(.Object, s), modelname, namespace)
		}
		return(.Object)
})

setMethod("genericExpRename", signature("MxExpectationBA81"),
	function(.Object, oldname, newname) {
		.Object@item <- renameReference(.Object@item, oldname, newname)
		.Object@mean <- renameReference(.Object@mean, oldname, newname)
		.Object@cov <- renameReference(.Object@cov, oldname, newname)
		.Object@data <- renameReference(.Object@data, oldname, newname)
		return(.Object)
})

##' Create a Bock & Aitkin (1981) expectation
##'
##' Used in conjuction with \link{mxFitFunctionML}, this expectation
##' models ordinal data with a modest number of latent dimensions.
##' Currently, only a multivariate Normal latent distribution is
##' supported.  An equal-interval quadrature is used to integrate over
##' the latent distribution.  When all items use the graded response
##' model and items are assumed conditionally independent then item
##' factor analysis is equivalent to a factor model.
##' 
##' The conditional likelihood of response \eqn{x_{ij}}{x[i,j]} to
##' item \eqn{j} from person \eqn{i} with item parameters
##' \eqn{\xi_j}{\xi[j]} and latent ability \eqn{\theta_i}{\theta[i]} is
##'
##' \deqn{L(x_i|\xi,\theta_i) = \prod_j \mathrm{Pr}(\mathrm{pick}=x_{ij} | \xi_j,\theta_i).}{%
##' L(x[i]|\xi,\theta[i]) = \prod_j Pr(pick=x[i,j] | \xi[j],\theta[i]).}
##'
##' Items are assumed conditionally independent.
##' That is, the outcome of one item is assumpted to not influence
##' another item after controlling for \eqn{\xi} and \eqn{\theta_i}{\theta[i]}.
##' The unconditional likelihood is obtained by integrating over
##' the latent distribution \eqn{\theta_i}{\theta[i]},
##'
##' \deqn{L(x_i|\xi) = \int L(x_i|\xi, \theta_i) L(\theta_i) \mathrm{d}\theta_i.}{%
##' L(x[i]|\xi) = \int L(x[i]|\xi, \theta[i]) L(\theta[i]) d \theta[i].}
##'
##' With an assumption that examinees are independently and identically distributed,
##' we can sum the individual log likelihoods,
##'
##' \deqn{\mathcal{L}=\sum_i \log L(x_i | \xi).}{%
##' L=\sum_i \log L(x[i] | \xi).}
##'
##' Response models \eqn{\mathrm{Pr}(\mathrm{pick}=x_{ij} |
##' \xi_j,\theta_i)}{Pr(pick=x[i,j] | \xi[j],\theta[i])}
##' are not implemented in OpenMx, but are imported
##' from the \href{https://cran.r-project.org/package=rpf}{RPF}
##' package. You must pass a list of models obtained from the RPF
##' package in the `ItemSpec' argument. All item models must use the
##' same number of latent factors although some of these factor
##' loadings can be constrained to zero in the item parameter matrix.
##' The `item' matrix contains item parameters with one item per
##' column in the same order at ItemSpec.
##'
##' The `qpoints' and `qwidth' argument control the fineness and
##' width, respectively, of the equal-interval quadrature grid.  The
##' integer `qpoints' is the number of points per dimension. The
##' quadrature extends from negative qwidth to positive qwidth for
##' each dimension. Since the latent distribution defaults to standard
##' Normal, qwidth can be regarded as a value in Z-score units.
##'
##' The optional `mean' and `cov' arguments permit modeling of the
##' latent distribution in multigroup models (in a single group, the
##' latent distribution must be fixed). A separate latent covariance
##' model is used in combination with mxExpectationBA81. The point
##' mass distribution contained in the quadrature is converted into a
##' multivariate Normal distribution by
##' \link{mxDataDynamic}. Typically \link{mxExpectationNormal} is used
##' to fit a multivariate Normal model to these data. Some intricate
##' programming is required.  Examples are given in the manual.
##' mxExpectationBA81 uses a sample size of \eqn{N} for the covariance
##' matrix. This differs from \link{mxExpectationNormal} which uses a
##' sample size of \eqn{N-1}.
##'
##' The `verbose' argument enables diagnostics that are mainly of
##' interest to developers.
##'
##' When a two-tier covariance matrix is recognized, this expectation
##' automatically enables analytic dimension reduction (Cai, 2010).
##' 
##' The optional `weightColumn' is superceded by the weight
##' argument in \link{mxData}. For data with many repeated
##' response patterns, model evaluation time can be
##' reduced. An easy way to transform your data into this form is to
##' use \link[rpf]{compressDataFrame}. Non-integer weights are supported except for
##' \link[rpf]{EAPscores}.
##'
##' mxExpectationBA81 requires \link{mxComputeEM}. During a typical
##' optimization run, latent abilities are assumed for examinees
##' during the E-step.  These examinee scores are implied by the
##' previous iteration's parameter vector. This can be overridden
##' using the `EstepItem' argument.  This is mainly of use to
##' developers for checking item parameter derivatives.
##'
##' Common univariate priors are available from
##' \link[ifaTools]{univariatePrior}.  The standard Normal
##' distribution of the quadrature acts like a prior distribution for
##' difficulty. It is not necessary to impose any additional Bayesian
##' prior on difficulty estimates (Baker & Kim, 2004, p. 196).
##'
##' Many estimators are available for standard errors. Oakes is
##' recommended (see \link{mxComputeEM}).  Also available are
##' Supplement EM (\link{mxComputeEM}), Richardson extrapolation
##' (\link{mxComputeNumericDeriv}), likelihood-based confidence
##' intervals (\link{mxCI}), and the covariance of the rowwise
##' gradients.
##'
##' @aliases MxExpectationBA81-class show,MxExpectationBA81-method print,MxExpectationBA81-method
##' @param ItemSpec a single item model (to replicate) or a list of
##' item models in the same order as the column of \code{ItemParam}
##' @param item the name of the mxMatrix holding item parameters
##' with one column for each item model with parameters starting at
##' row 1 and extra rows filled with NA
##' @param ...  Not used.  Forces remaining arguments to be specified by name.
##' @param qpoints number of points to use for equal interval quadrature integration (default 49L)
##' @param qwidth the width of the quadrature as a positive Z score (default 6.0)
##' @param mean the name of the mxMatrix holding the mean vector
##' @param cov the name of the mxMatrix holding the covariance matrix
##' @param verbose the level of runtime diagnostics (default 0L)
##' @param weightColumn the name of the column in the data containing the row weights (DEPRECATED)
##' @param EstepItem a simple matrix of item parameters for the
##' E-step. This option is mainly of use for debugging derivatives.
##' @param debugInternal when enabled, some of the internal tables are
##' returned in $debug. This is mainly of use to developers.
##' @seealso \href{https://cran.r-project.org/package=rpf}{RPF}
##' @references
##' Bock, R. D., & Aitkin, M. (1981). Marginal maximum likelihood estimation of item
##' parameters: Application of an EM algorithm. \emph{Psychometrika, 46}, 443-459.
##'
##' Cai, L. (2010). A two-tier full-information item factor analysis
##' model with applications. \emph{Psychometrika, 75}, 581-612.
##'
##' Pritikin, J. N., Hunter, M. D., & Boker, S. M. (2015). Modular
##' open-source software for Item Factor Analysis. \emph{Educational and
##' Psychological Measurement, 75}(3), 458-474
##'
##' Pritikin, J. N. & Schmidt, K. M. (in press). Model builder for
##' Item Factor Analysis with OpenMx. \emph{R Journal}.
##' 
##' Seong, T. J. (1990). Sensitivity of marginal maximum likelihood
##' estimation of item and ability parameters to the characteristics
##' of the prior ability distributions. \emph{Applied Psychological
##' Measurement, 14}(3), 299-311.
##'
##' @examples
##' library(OpenMx)
##' library(rpf)
##' 
##' numItems <- 14
##' 
##' # Create item specifications
##' spec <- list()
##' for (ix in 1:numItems) { spec[[ix]] <- rpf.grm(outcomes=sample(2:7, 1)) }
##' names(spec) <- paste("i", 1:numItems, sep="")
##' 
##' # Generate some random "true" parameter values
##' correct.mat <- mxSimplify2Array(lapply(spec, rpf.rparam))
##' 
##' # Generate some example data
##' data <- rpf.sample(500, spec, correct.mat)
##' 
##' # Create a matrix of item parameters with starting values
##' imat <- mxMatrix(name="item",
##'                  values=mxSimplify2Array(lapply(spec, rpf.rparam)))
##' rownames(imat)[1] <- 'f1'
##' imat$free[!is.na(correct.mat)] <- TRUE
##' imat$values[!imat$free] <- NA
##' 
##' # Create a compute plan
##' plan <- mxComputeSequence(list(
##'   mxComputeEM('expectation', 'scores',
##'               mxComputeNewtonRaphson(), information="oakes1999",
##'               infoArgs=list(fitfunction='fitfunction')),
##'   mxComputeHessianQuality(),
##'   mxComputeStandardError(),
##'   mxComputeReportDeriv()))
##' 
##' # Build the OpenMx model
##' grmModel <- mxModel(model="grm1", imat,
##'                     mxData(observed=data, type="raw"),
##'                     mxExpectationBA81(ItemSpec=spec),
##'                     mxFitFunctionML(),
##'                     plan)
##' 
##' grmModel <- mxRun(grmModel)
##' summary(grmModel)
mxExpectationBA81 <- function(ItemSpec, item="item", ...,
			      qpoints=49L, qwidth=6.0, mean="mean", cov="cov",
			      verbose=0L, weightColumn=NA_integer_,
			      EstepItem=NULL, debugInternal=FALSE) {

	if (length(list(...)) > 0) {
		stop(paste("Remaining parameters must be passed by name", deparse(list(...))))
	}

	if (packageVersion("rpf") < "0.28") stop("Please install 'rpf' version 0.28 or newer")
	if (qpoints < 3) {
		stop("qpoints must be 3 or greater")
	}
	if (missing(qwidth)) qwidth <- 6
	if (qwidth <= 0) {
		stop("qwidth must be positive")
	}
  
	if (!is.list(ItemSpec)) ItemSpec <- list(ItemSpec)

	minItemsPerScore <- 0L # remove TODO

	if (is.na(weightColumn)) weightColumn <- as.integer(weightColumn)

	return(new("MxExpectationBA81", ItemSpec, item, EstepItem,
		   qpoints, qwidth, mean, cov, as.integer(verbose), debugInternal,
		   minItemsPerScore, weightColumn))
}

##' Like simplify2array but works with vectors of different lengths
##'
##' Vectors are filled column-by-column into a matrix. Shorter vectors
##' are padded with NAs to fill whole columns.
##' @param x a list of vectors
##' @param higher whether to produce a higher rank array (defaults to FALSE)
##' @examples
##' v1 <- 1:3
##' v2 <- 4:5
##' v3 <- 6:10
##' mxSimplify2Array(list(v1,v2,v3))
##'
##' #     [,1] [,2] [,3]
##' # [1,]    1    4    6
##' # [2,]    2    5    7
##' # [3,]    3   NA    8
##' # [4,]   NA   NA    9
##' # [5,]   NA   NA   10

mxSimplify2Array <- function(x, higher=FALSE) {
	if (higher) {
		stop("higher=TRUE is not implemented. Consider using simplify2array")
	}
  len <- sapply(x, length)
  biggest <- which(len == max(len))[1]
  out <- matrix(NA, nrow=max(len), ncol=length(x))
  for (iter in 1:length(x)) {
	  if (len[iter] == 0) next
    out[1:len[iter],iter] <- x[[iter]]
  }
  colnames(out) <- names(x)
  rownames(out) <- names(x[[biggest]])
  out
}
