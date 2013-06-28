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
	   ItemSpec = "MxCharOrNumber",
	   Design = "MxOptionalCharOrNumber",
           ItemParam = "MxCharOrNumber",
           EItemParam = "MxOptionalCharOrNumber",
           CustomPrior = "MxOptionalCharOrNumber",
	   qpoints = "numeric",
	   qwidth = "numeric",
	   cache = "logical",
	   rescale = "logical",
	   scores = "character",
	   mean = "MxCharOrNumber",
	   cov = "MxCharOrNumber"),
         contains = "MxBaseExpectation")

setMethod("initialize", "MxExpectationBA81",
          function(.Object, ItemSpec, ItemParam, EItemParam, CustomPrior, Design,
		   qpoints, qwidth, cache, mean, cov, rescale, scores,
		   name = 'expectation') {
            .Object@name <- name
            .Object@ItemSpec <- ItemSpec
            .Object@Design <- Design
            .Object@ItemParam <- ItemParam
            .Object@EItemParam <- EItemParam
            .Object@qpoints <- qpoints
            .Object@qwidth <- qwidth
            .Object@cache <- cache
            .Object@rescale <- rescale
            .Object@scores <- scores
            .Object@CustomPrior <- CustomPrior
            .Object@data <- as.integer(NA)
	    .Object@mean <- mean
	    .Object@cov <- cov
            return(.Object)
          }
)

setMethod("genericExpDependencies", signature("MxExpectationBA81"),
	  function(.Object, dependencies) {
		  sources <- c(.Object@CustomPrior, .Object@ItemParam,
			       .Object@EItemParam,
			       .Object@mean, .Object@cov,
			       .Object@ItemSpec, .Object@Design)
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
		  for (s in c("data", "ItemParam", "EItemParam", "CustomPrior",
			      "ItemSpec", "Design", "mean", "cov")) {
			  if (is.null(slot(.Object, s))) next;
			  slot(.Object, s) <-
			    imxLocateIndex(flatModel, slot(.Object, s), name)
		  }

					# How to get the data object?
		  ## if (.Object@data@type != 'raw') {
		  ##   stop(paste(typeof(.Object), "only supports raw data"));
		  ## }
		  return(.Object)
	  })

setMethod("genericExpFunNamespace", signature("MxExpectationBA81"), 
	function(.Object, modelname, namespace) {
		.Object@name <- imxIdentifier(modelname, .Object@name)
		for (s in c("ItemParam", "EItemParam", "CustomPrior", "ItemSpec",
			    "Design", "mean", "cov")) {
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
##' @param ItemSpec one column for each item containing the model ID
##' (see \code{\link{mxLookupIRTItemModelID}}), numDimensions, and
##' numOutcomes (3 rows)
##' @param ItemParam one column for each item with parameters starting
##' at row 1 and extra rows filled with NA
##' @param Design one column per item assignment of person abilities
##' to item dimensions (optional)
##' @param qpoints number of points to use for rectangular quadrature integrations (default 49)
##' See Seong (1990) for some considerations on specifying this parameter.
##' @param cache whether to cache part of the expectation calculation
##' (enabled by default).
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

mxExpectationBA81 <- function(ItemSpec, ItemParam, EItemParam=NULL, CustomPrior=NULL, Design=NULL,
			      qpoints=NULL, qwidth=6.0, cache=TRUE, mean=NULL, cov=NULL,
			      rescale=TRUE, scores="omit") {

	if (missing(qpoints)) qpoints <- 49
	if (qpoints < 3) {
		stop("qpoints should be 3 or greater")
	}
	if (missing(qwidth)) qwidth <- 6
	if (qwidth <= 0) {
		stop("qwidth must be positive")
	}
  
	score.options <- c("omit", "unique", "full")
	if (!match(scores, score.options)) {
		stop(paste("Valid score options are", deparse(score.options)))
	}

	return(new("MxExpectationBA81", ItemSpec, ItemParam, EItemParam, CustomPrior, Design,
		   qpoints, qwidth, cache, mean, cov, rescale, scores))
}
