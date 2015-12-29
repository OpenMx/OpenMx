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
	     dataColumns="numeric",
	   dims = "character",
	   verbose = "integer",
	     minItemsPerScore = "integer",
	     weightColumn = "MxCharOrNumber"),
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

		  .Object@dataColumns <- as.numeric(dc)
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
          # not sure what goes here yet
		return(.Object)
})

##' Create a Bock & Aitkin (1981) expectation
##'
##' When a two-tier covariance matrix is recognized, this expectation
##' automatically enables analytic dimension reduction (Cai, 2010).
##' 
##' The standard Normal distribution of the quadrature acts like a
##' prior distribution for difficulty. It is not necessary to impose
##' any additional Bayesian prior on difficulty estimates (Baker &
##' Kim, 2004, p. 196).
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
##' @param weightColumn the name of the column in the data containing the row weights (default NA)
##' @param EstepItem a simple matrix of item parameters for the
##' E-step. This option is mainly of use for debugging derivatives.
##' @param debugInternal when enabled, some of the internal tables are
##' returned in $debug. This is mainly of use to developers.
##' @seealso \href{http://cran.r-project.org/package=rpf}{RPF}
##' @references
##' Bock, R. D., & Aitkin, M. (1981). Marginal maximum likelihood estimation of item
##' parameters: Application of an EM algorithm. Psychometrika, 46, 443-459.
##'
##' Cai, L. (2010). A two-tier full-information item factor analysis
##' model with applications. Psychometrika, 75, 581-612.
##'
##' Seong, T. J. (1990). Sensitivity of marginal maximum likelihood
##' estimation of item and ability parameters to the characteristics
##' of the prior ability distributions. Applied Psychological
##' Measurement, 14(3), 299-311.

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
