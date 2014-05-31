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
	   ItemParam = "MxCharOrNumber",
	   EItemParam = "MxOptionalMatrix",
	   qpoints = "numeric",
	   qwidth = "numeric",
	   scores = "character",
	   mean = "MxCharOrNumber",
	   cov = "MxCharOrNumber",
	     debugInternal="logical",
	     dataColumns="integer",
	   dims = "character",
	   numStats = "numeric",
	   verbose = "integer",
	     output = "list",
	     debug = "list",
	     naAction = "character",
	     minItemsPerScore = "integer",
	     weightColumn = "MxCharOrNumber"),
         contains = "MxBaseExpectation")

setMethod("initialize", "MxExpectationBA81",
          function(.Object, ItemSpec, ItemParam, EItemParam,
		   qpoints, qwidth, mean, cov, scores, verbose, debugInternal,
		   naAction, minItemsPerScore, weightColumn, name = 'expectation') {
            .Object@name <- name
	    .Object@ItemSpec <- ItemSpec
	    .Object@ItemParam <- ItemParam
	    .Object@EItemParam <- EItemParam
            .Object@qpoints <- qpoints
            .Object@qwidth <- qwidth
            .Object@scores <- scores
            .Object@data <- as.integer(NA)
	    .Object@mean <- mean
	    .Object@cov <- cov
	    .Object@verbose <- verbose
	    .Object@debugInternal <- debugInternal
	    .Object@naAction <- naAction
	    .Object@minItemsPerScore <- minItemsPerScore
	    .Object@weightColumn <- weightColumn
            return(.Object)
          }
)

setMethod("genericExpDependencies", signature("MxExpectationBA81"),
	  function(.Object, dependencies) {
		  sources <- c(.Object@mean, .Object@cov,
			       .Object@ItemParam)
		  dependencies <- imxAddDependency(sources, .Object@name, dependencies)
		  return(dependencies)
	  })

setMethod("genericExpFunConvert", signature("MxExpectationBA81"), 
	  function(.Object, flatModel, model, labelsData, defVars, dependencies) {
		  modelname <- imxReverseIdentifier(model, .Object@name)[[1]]
		  if(is.na(.Object@data)) {
			  msg <- paste(typeof(.Object),
				       "does not have a dataset associated with it in model",
				       omxQuotes(modelname))
			  stop(msg, call.=FALSE)
		  }

		  verifyMvnNames(.Object@cov, .Object@mean, "prior", flatModel, model@name, class(.Object))

		  fnames <- dimnames(flatModel[[.Object@cov]])[[1]]

		  ItemParam <- flatModel@matrices[[.Object@ItemParam]]
		  if (is.null(rownames(ItemParam)) || any(rownames(ItemParam)[1:length(fnames)] != fnames)) {
			  msg <- paste("The first", length(fnames), "rownames of ItemParam",
				       "must be", omxQuotes(fnames))
			  stop(msg)
		  }

		  name <- .Object@name
		  for (s in c("data", "ItemParam", "mean", "cov")) {
			  if (is.null(slot(.Object, s))) next
			  slot(.Object, s) <-
			    imxLocateIndex(flatModel, slot(.Object, s), name)
		  }

		  .Object@dims <- colnames(flatModel@datasets[[.Object@data + 1]]@observed) # for summary

		  if (is.null(colnames(ItemParam))) {
			  stop(paste(class(.Object),
				     ": the ItemParam matrix must have column names",
				     "to match against the observed data column names."))
		  } else {
			  dc <- match(colnames(ItemParam), .Object@dims) - 1L
			  if (any(is.na(dc))) {
				  msg <- paste("Some items are not found among the observed data columns:",
					       omxQuotes(colnames(ItemParam)[is.na(dc)]))
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
				  .Object@weightColumn <- wc
			  }
		  }
		  return(.Object)
	  })

setMethod("qualifyNames", signature("MxExpectationBA81"), 
	function(.Object, modelname, namespace) {
		.Object@name <- imxIdentifier(modelname, .Object@name)
		for (s in c("ItemParam", "mean", "cov")) {
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
##' The standard Normal distribution of the quadrature acts like a
##' prior distribution for difficulty. It is not necessary to impose
##' any additional Bayesian prior on difficulty estimates (Baker &
##' Kim, 2004, p. 196).
##' 
##' @param ItemParam one column for each item with parameters starting
##' at row 1 and extra rows filled with NA
##' @param qpoints number of points to use for rectangular quadrature integrations (default 49)
##' See Seong (1990) for some considerations on specifying this parameter.
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

mxExpectationBA81 <- function(ItemSpec, ItemParam,
			      qpoints=NULL, qwidth=6.0, mean=NULL, cov=NULL,
			      scores="omit", verbose=0L, EItemParam=NULL, debugInternal=FALSE,
			      naAction="fail", minItemsPerScore=1L, weightColumn=NA_integer_) {

	if (packageVersion("rpf") < "0.28") stop("Please install 'rpf' version 0.28 or newer")
	if (missing(qpoints)) qpoints <- 49
	if (qpoints < 3) {
		stop("qpoints must be 3 or greater")
	}
	if (missing(qwidth)) qwidth <- 6
	if (qwidth <= 0) {
		stop("qwidth must be positive")
	}
  
	score.options <- c("omit", "full")
	if (!match(scores, score.options)) {
		stop(paste("Valid score options are", omxQuotes(score.options)))
	}

	if (!is.list(ItemSpec)) ItemSpec <- list(ItemSpec)

	actions <- c("pass", "fail")
	if (!match(naAction, actions)) stop(paste("Valid choices for naAction are", omxQuotes(actions)))

	if (minItemsPerScore < 0L) stop("minItemsPerScore must be non-negative")

	if (is.na(weightColumn)) weightColumn <- as.integer(weightColumn)

	return(new("MxExpectationBA81", ItemSpec, ItemParam, EItemParam,
		   qpoints, qwidth, mean, cov, scores, verbose, debugInternal,
		   naAction, minItemsPerScore, weightColumn))
}
