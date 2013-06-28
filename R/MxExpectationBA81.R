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
           CustomPrior = "MxOptionalCharOrNumber",
	   GHpoints = "numeric",  # rename, not necessarily G-H TODO
	   GHarea = "numeric",
	   cache = "logical",
	   doRescale = "logical"),
         contains = "MxBaseExpectation")

setMethod("initialize", "MxExpectationBA81",
          function(.Object, ItemSpec, ItemParam, CustomPrior, Design, RPF,
		   quadrature, cache, doRescale, name = 'expectation') {
            .Object@name <- name
            .Object@ItemSpec <- ItemSpec
            .Object@Design <- Design
            .Object@ItemParam <- ItemParam
            .Object@RPF <- RPF
            .Object@GHpoints <- quadrature[[1]]
            .Object@GHarea <- quadrature[[2]]
            .Object@cache <- cache
            .Object@CustomPrior <- CustomPrior
            .Object@data <- as.integer(NA)
	    .Object@doRescale <- doRescale
            return(.Object)
          }
)

setMethod("genericExpDependencies", signature("MxExpectationBA81"),
	  function(.Object, dependencies) {
		  sources <- c(.Object@CustomPrior, .Object@ItemParam,
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
		  for (s in c("data", "ItemParam", "CustomPrior", "ItemSpec", "Design")) {
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
		for (s in c("ItemParam", "CustomPrior", "ItemSpec", "Design")) {
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
##' specific dimensions.
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

mxExpectationBA81 <- function(ItemSpec, ItemParam, CustomPrior=NULL, Design=NULL,
			      RPF=NULL, GHpoints=NULL, cache=TRUE, quadrature=NULL, doRescale=TRUE) {
	if (missing(quadrature)) {
		if (missing(GHpoints)) GHpoints <- 10
		if (GHpoints < 3) {
			error("GHpoints should be 3 or greater")
		}
		quadrature <- rpf.GaussHermiteData(GHpoints)
	}
  
	return(new("MxExpectationBA81", ItemSpec, ItemParam, CustomPrior, Design, RPF,
		   quadrature, cache, doRescale))
}

imxUniformQuadratureData <- function(n, width) {
  x <- seq(-width, width, length.out=n)
  list(x=x, w=rep(1/n, n))
}

imxEqualIntervalQuadratureData <- function(n, width) {
  x <- seq(-width, width, length.out=n)
  list(x=x, w=dnorm(x)/sum(dnorm(x)))
}
