#
#   Copyright 2007-2016 The OpenMx Project
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
	    list("Default optimizer" = imxDetermineDefaultOptimizer())),
	'mxByrow' = FALSE,
	'mxCondenseMatrixSlots' = FALSE,
	'mxShowDimnames' = TRUE,
	'mxPrintUnitTests' = TRUE,
	'swift.initialexpr' = "library(OpenMx)")
}

##' imxHasOpenMP
##'
##' This is an internal function exported for those people who know
##' what they are doing.
imxHasOpenMP <- function() .Call(hasOpenMP_wrapper)

.onLoad <- function(libname, pkgname) {
	# http://stackoverflow.com/questions/12598242/global-variables-in-packages-in-r
	assign("pkg_globals", new.env(), envir=parent.env(environment()))
	mxSetDefaultOptions()
}

.onAttach <- function(libname, pkgname) {
	if (!imxHasOpenMP()) {
		packageStartupMessage("OpenMx is not compiled to take advantage of computers with multiple cores.")
	}
	if (!is.na(match("package:expm", search()))) {
		packageStartupMessage(paste("** Holy cannoli! You must be a pretty advanced and awesome user.",
					    "The expm package is loaded.",
					    "Note that expm defines %^% as repeated matrix multiplication (matrix to a power)",
					    "whereas OpenMx defines the same operation as",
					    "elementwise powering of one matrix by another (Kronecker power)."))
	}
}

#' Test thread-safe output code
#'
#' This is the code that the backend uses to write diagnostic
#' information to standard error. This function should not be called
#' from R. We make it available only for testing.
#' 
#' @param str the character string to output
imxLog <- function(str) .Call(Log_wrapper, str)

#' OpenMx: An package for Structural Equation Modeling and Matrix Algebra Optimization
#'
#' OpenMx is a package for structural equation modeling, matrix algebra optimization and other statistical estimation problems. 
#' Try the example below. We try and have useful help files: for instance help(\code{\link{mxRun}}) to learn more. Also the reference manual
#'
#' @details OpenMx solves algebra optimization and statistical estimation problems using matrix algebra. 
#' Most users use it for Structural equation modeling.
#' 
#' The core function is \code{\link{mxModel}}, which makes a model. Models are containers for data, matrices, \code{\link{mxPath}}s
#' algebras, bounds and constraints. Models most often have an expectation function (e.g., \code{\link{mxExpectationNormal}}) 
#' to calculate the expectations for the model. Models need a fit function. Several of these are built-in (e.g., \link{mxFitFunctionML} )
#' OpenMx also allows user-defined fit functions for purposes not covered by the built-in functions. (e.g., \link{mxFitFunctionR} or \link{mxFitFunctionAlgebra}).
#' 
#' Once built, the resulting mxModel can be run (i.e., optimized) using  \code{\link{mxRun}}. This returns the fitted model.
#' 
#' You can see the resulting parameter estimates, algebra evaluation etc using \code{\link[=summary.MxModel]{summary}}(yourModel)
#' 
#' The user's manual is online (see reference below), but functions \link{mxRun}, \link{mxModel}, \link{mxMatrix}
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
#' Michael C. Neale, Michael D. Hunter, Joshua N. Pritikin, Mahsa Zahery, Timothy R. Brick Robert M.
#' Kickpatrick, Ryne Estabrook, Timothy C. Bates, Hermine H. Maes, Steven M. Boker. (in press).
#' OpenMx 2.0: Extended structural equation and statistical modeling. \emph{Psychometrika}. 
#' DOI: 10.1007/s11336-014-9435-8
#' 
#' Steven M. Boker, Michael C. Neale, Hermine H. Maes, Michael J. Wilde, Michael Spiegel, Timothy R. Brick,
#' Jeffrey Spies, Ryne Estabrook, Sarah Kenny, Timothy C. Bates, Paras Mehta, and John Fox. (2011)
#' OpenMx: An Open Source Extended Structural Equation Modeling Framework. 
#' \emph{Psychometrika}, 306-317. DOI:10.1007/s11336-010-9200-6
#' 
#' Steven M. Boker, Michael C. Neale, Hermine H. Maes, Michael J. Wilde, Michael Spiegel, Timothy R. Brick, Ryne
#' Estabrook, Timothy C. Bates, Paras Mehta, Timo von Oertzen, Ross J. Gore, Michael D. Hunter, Daniel C.
#' Hackett, Julian Karch, Andreas M. Brandmaier, Joshua N. Pritikin, Mahsa Zahery, Robert M. Kirkpatrick, 
#' Yang Wang, and Charles Driver. (2016) OpenMx 2 User Guide. 
#' http://openmx.psyc.virginia.edu/docs/OpenMx/latest/OpenMxUserGuide.pdf
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
