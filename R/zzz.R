#
#   Copyright 2007-2014 The OpenMx Project
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

##' Reset global options to the default
mxSetDefaultOptions <- function() {
   options('mxDefaultType' = 'default', 
	'mxOptions' = c(npsolOptions, checkpointOptions, otherOptions,
	    list("Default optimizer" = determineDefaultOptimizer())),
	'mxByrow' = FALSE,
	'mxShowDimnames' = TRUE,
	'mxPrintUnitTests' = TRUE,
	'swift.initialexpr' = "library(OpenMx)")
}

.onLoad <- function(libname, pkgname) mxSetDefaultOptions()

#' OpenMx: An package for Structural Equation Modeling and Matrix Algebra Optimization
#'
#' OpenMx is a package for structural equation modeling, matrix algebra optimization and other statistical estimation problems. 
#' Try the example below. Have a look at the references and at functions like \code{\link{umxRun}} to learn more.
#'
#' @details OpenMx solves algebra optimization and statistical estimation problems using matrix algebra. 
#' Most users use it for Structural equation modeling.
#' 
#' The core function is \code{\link{mxModel}}, which makes a model. Models are containers for data, matrices, \code{\link{mxPath}}s
#' algebras, bounds and constraints. Models most often have an expectation function (e.g., \code{\link{mxExpectationNormal}}) 
#' to calculate the expectations for the model. Models need a fit function. Several of these are built-in (e.g., \link{mxFitFunctionML} )
#' OpenMx also allows user-defined fit functions for purposes not covered by the built-in functions. (e.g., \link{mxFitFunctionR} or \link{mxFitFunctionAlgebra}).
#' 
#' Once built, the resulting mxModel can be run (i.e., optimized) using \link{mxRun}. This returns the fitted model.
#' 
#' You can see the resulting parameter estimates, algebra evaluation etc using summary(yourModel)
#' 
#' The user\u2019s manual is online (see reference below), but functions \link{mxRun}, \link{mxModel}, \link{mxMatrix}
#' all have working examples to get you started as well.
#'
#' The main OpenMx functions are: \link{mxAlgebra}, \link{mxBounds}, \link{mxCI}, \link{mxConstraint}, \link{mxData}, 
#' \link{mxMatrix}, \link{mxModel}, and \link{mxPath}
#' 
#' Expectation functions include \link{mxExpectationNormal}, \link{mxExpectationRAM}, \link{mxExpectationLISREL}, and \link{mxExpectationStateSpace};
#' 
#' Fit functions include \link{mxFitFunctionML}, \link{mxFitFunctionAlgebra}, \link{mxFitFunctionRow} and \link{mxFitFunctionR}.
#' 
#' OpenMx comes with several useful datasets built-in. Access them using data(package="OpenMx")
#' 
#' 
#' To cite package 'OpenMx' in publications use:
#' 
#' Steven M. Boker, Michael C. Neale, Hermine H. Maes, Michael J. Wilde, Michael Spiegel, Timothy R. Brick,
#' Jeffrey Spies, Ryne Estabrook, Sarah Kenny, Timothy C. Bates, Paras Mehta, and John Fox. (2011) OpenMx: An
#' Open Source Extended Structural Equation Modeling Framework. Psychometrika.
#' 
#' Steven M. Boker, Michael C. Neale, Hermine H. Maes, Michael J. Wilde, Michael Spiegel, Timothy R. Brick, Ryne
#' Estabrook, Timothy C. Bates, Paras Mehta, Timo von Oertzen, Ross J. Gore, Michael D. Hunter, Daniel C.
#' Hackett, Julian Karch and Andreas M. Brandmaier. (2014) OpenMx 2 User Guide.
#' 
#' @references The OpenMx User's guide can be found at \url{http://openmx.psyc.virginia.edu/documentation}
#'
#' @examples
#' library(OpenMx)
#' data(demoOneFactor)
#' # ===============================
#' # = Make and run a 1-factor CFA =
#' # ===============================
#' 
#' latents  = c("G") # the latent factor
#' manifests = names(demoOneFactor) # manifest variables to be modeled
#' # ====================
#' # = Make the MxModel =
#' # ====================
#' m1 <- mxModel("One Factor", type = "RAM", 
#' 	manifestVars = manifests, latentVars = latents, 
#' 	mxPath(from = latents, to = manifests),
#' 	mxPath(from = manifests, arrows = 2),
#' 	mxPath(from = latents, arrows = 2, free = FALSE, values = 1.0),
#' 	mxData(cov(demoOneFactor), type = "cov", numObs = 500)
#' )
#' 
#' # ===============================
#' # = mxRun it and get a summary! =
#' # ===============================
#' 
#' m1 = mxRun(m1)
#' summary(m1, show = "std")
#'
#' @docType package
#' @name OpenMx
NULL
