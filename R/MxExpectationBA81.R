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


setClassUnion("MxOptionalFunction", c("NULL", "function"))

setClass(Class = "MxExpectationBA81",
         representation = representation(
	   ItemSpec = "MxCharOrNumber",
	   Design = "MxOptionalCharOrNumber",
           ItemParam = "MxCharOrNumber",
	   RPF = "MxOptionalFunction",
           ItemPrior = "MxCharOrNumber",
	   GHpoints = "numeric",
	   GHarea = "numeric",
	   cache = "logical"),
         contains = "MxBaseExpectation")

setMethod("initialize", "MxExpectationBA81",
          function(.Object, ItemSpec, ItemParam, ItemPrior, Design, RPF, gh, cache, name = 'expectation') {
            .Object@name <- name
            .Object@ItemSpec <- ItemSpec
            .Object@Design <- Design
            .Object@ItemParam <- ItemParam
            .Object@RPF <- RPF
            .Object@GHpoints <- gh[[1]]
            .Object@GHarea <- gh[[2]]
            .Object@cache <- cache
            .Object@ItemPrior <- ItemPrior
            .Object@data <- as.integer(NA)
            return(.Object)
          }
)

setMethod("genericExpDependencies", signature("MxExpectationBA81"),
	  function(.Object, dependencies) {
		  sources <- c(.Object@ItemPrior, .Object@ItemParam,
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
		  for (s in c("data", "ItemParam", "ItemPrior", "ItemSpec", "Design")) {
			  if (is.null(slot(.Object, s))) next;
			  slot(.Object, s) <-
			    imxLocateIndex(flatModel, slot(.Object, s), name)
		  }

					# How to get the data object?
		  ## if (.Object@data@type != 'raw') {
		  ##   error(paste(typeof(.Object), "only supports raw data"));
		  ## }
		  return(.Object)
	  })

setMethod("genericExpFunNamespace", signature("MxExpectationBA81"), 
	function(.Object, modelname, namespace) {
		.Object@name <- imxIdentifier(modelname, .Object@name)
		for (s in c("ItemParam", "ItemPrior", "ItemSpec", "Design")) {
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
##' The RPF function is an optional argument for experimenting with
##' novel response probability functions. The model ID in the ItemSpec
##' is ignored if you provide your own function. There is no provision
##' to provide analytic gradients from R. If you want to provide
##' analytic gradients then you need to implement them in C.
##'
##' The standard Normal distribution of the quadrature acts like a
##' prior distribution for difficulty. It is not necessary to impose
##' any additional Bayesian prior on difficulty estimates (Baker &
##' Kim, 2004, p. 196).
##' 
##' The cache more than halves CPU time, especially if you have many
##' specific dimensions. With the cache enabled, memory usage is about
##' 8 * (numItems + numUnique * (2 + numThreads * (numSpecific+1)) +
##' numUnique * product(quadrature grids) * numSpecific). With the
##' cache disabled, memory use is limited to about 8 * (numItems + 2 *
##' numUnique + (numUnique+1) * (numThreads * (numSpecific+1)) bytes.
##' For example, suppose your model has 20 items, 500 unique response
##' patterns, 2 primary dimensions, 3 specific dimensions, a 40 point
##' quadrature and your computer offers 8 way multiprocessing. With
##' the cache enabled, 8 * (20 + 500 * (2+8 * (3+1)) + 500 * 40^3 * 3)
##' = 733 MiB. With the cache disabled, 8 * (20 + 2 * 500 + (500 + 1) *
##' (8 * (3 + 1))) = 134 KiB.
##'
##' @param ItemSpec one column for each item containing the model ID
##' (see \code{\link{mxLookupIRTItemModelID}}), numDimensions, and
##' numOutcomes (3 rows)
##' @param ItemParam one column for each item with parameters starting
##' at row 1 and extra rows filled with NA
##' @param Design one column per item assignment of person abilities
##' to item dimensions (optional)
##' @param RPF a function to turn item parameters into a response
##' probability matrix (optional)
##' @param GHpoints number of points to use for Gauss-Hermite quadrature integrations (default 10)
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

mxExpectationBA81 <- function(ItemSpec, ItemParam, ItemPrior, Design=NULL, RPF=NULL, GHpoints=10, cache=TRUE) {
  if (GHpoints < 3) {
    error("GHpoints should be 3 or greater")
  }
  gh <- imxGaussHermiteData(GHpoints)
  
  return(new("MxExpectationBA81", ItemSpec, ItemParam, ItemPrior, Design, RPF, gh, cache))
}

# Not used, just informational
imxAdhocGaussHermiteData <- function(n, width) {
  x <- seq(-width, width, length.out=n)
  list(x=x, w=dnorm(x)/sum(dnorm(x)))
}

imxGaussHermiteData <- function(n) .Call(omxGaussHermiteData, n)

##' Convert an IRT item model name to an ID
##'
##' drm1 is the standard 3PL. drm is the multidimensional version of
##' the 3PL. gpcm1 is Generalized Partial Credit Model.
##'
##' @param name name of the item model (string)
##' @return the integer ID assigned to the given model
mxLookupIRTItemModelID <- function(name) {
	models <- .Call(omx_get_rpf_names)
	match(name, models) - 1;
}
