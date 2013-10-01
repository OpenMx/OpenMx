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
	   CustomPrior = "MxOptionalCharOrNumber",
	   design = "MxOptionalMatrix",
	   qpoints = "numeric",
	   qwidth = "numeric",
	   scores = "character",
	   scores.out = "matrix",
	   mean = "MxCharOrNumber",
	   cov = "MxCharOrNumber",
	   empirical.mean = "numeric",
	   empirical.cov = "matrix",
	   patternLikelihood = "numeric",
	   em.expected = "numeric",
	   dims = "character",
	   numStats = "numeric",
	   verbose = "logical"),
         contains = "MxBaseExpectation")

setMethod("initialize", "MxExpectationBA81",
          function(.Object, ItemSpec, ItemParam, EItemParam, CustomPrior, design,
		   qpoints, qwidth, mean, cov, scores, verbose,
		   name = 'expectation') {
            .Object@name <- name
	    .Object@ItemSpec <- ItemSpec
	    .Object@ItemParam <- ItemParam
	    .Object@EItemParam <- EItemParam
	    .Object@CustomPrior <- CustomPrior
            .Object@design <- design
            .Object@qpoints <- qpoints
            .Object@qwidth <- qwidth
            .Object@scores <- scores
            .Object@data <- as.integer(NA)
	    .Object@mean <- mean
	    .Object@cov <- cov
	    .Object@scores.out <- matrix()
	    .Object@verbose <- verbose
            return(.Object)
          }
)

setMethod("genericExpDependencies", signature("MxExpectationBA81"),
	  function(.Object, dependencies) {
		  sources <- c(.Object@mean, .Object@cov,
			       .Object@ItemParam, .Object@CustomPrior)
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
		  name <- .Object@name
		  for (s in c("data", "ItemParam", "CustomPrior",
			      "mean", "cov")) {
			  if (is.null(slot(.Object, s))) next;
			  slot(.Object, s) <-
			    imxLocateIndex(flatModel, slot(.Object, s), name)
		  }
		  .Object@dims <- colnames(flatModel@datasets[[.Object@data + 1]]@observed)
		  return(.Object)
	  })

setMethod("qualifyNames", signature("MxExpectationBA81"), 
	function(.Object, modelname, namespace) {
		.Object@name <- imxIdentifier(modelname, .Object@name)
		for (s in c("ItemParam", "CustomPrior", "mean", "cov")) {
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
##' @param design one column per item, assignment of person abilities
##' to item dimensions (optional)
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

mxExpectationBA81 <- function(ItemSpec, ItemParam, CustomPrior=NULL, design=NULL,
			      qpoints=NULL, qwidth=6.0, mean=NULL, cov=NULL,
			      scores="omit", verbose=FALSE, EItemParam=NULL) {

	if (packageVersion("rpf") < "0.15") stop("Please install 'rpf' version 0.15 or newer")
	if (missing(qpoints)) qpoints <- 49
	if (qpoints < 3) {
		stop("qpoints should be 3 or greater")
	}
	if (qpoints %% 2 == 0) {
		warning(paste("An even number of qpoints can obtain a better than true fit",
			      "in a single group model; Pick an odd number of qpoints"))
	}
	if (missing(qwidth)) qwidth <- 6
	if (qwidth <= 0) {
		stop("qwidth must be positive")
	}
  
	score.options <- c("omit", "unique", "full")
	if (!match(scores, score.options)) {
		stop(paste("Valid score options are", deparse(score.options)))
	}

	if (!missing(design) && !is.integer(design)) {
		stop("Design must be an integer matrix")
	}

	return(new("MxExpectationBA81", ItemSpec, ItemParam, EItemParam, CustomPrior, design,
		   qpoints, qwidth, mean, cov, scores, verbose))
}
