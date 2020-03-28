#
#   Copyright 2007-2019 by the individuals mentioned in the source code history
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
  options(
    'mxDefaultType' = 'default',
    'mxOptions' = c(npsolOptions, checkpointOptions, otherOptions,
                    list("Default optimizer" = imxDetermineDefaultOptimizer()),
                    "Number of Threads" = imxGetNumThreads()),
    'mxByrow' = FALSE,
    'mxCondenseMatrixSlots' = FALSE,
    'swift.initialexpr' = "library(OpenMx)")
}

##' imxHasOpenMP
##'
##' This is an internal function exported for those people who know
##' what they are doing.
imxHasOpenMP <- function() .Call(hasOpenMP_wrapper)

# Don't use .onAttach except for packageStartupMessage,
# see https://github.com/OpenMx/OpenMx/issues/98
.onLoad <- function(libname, pkgname) {
  .Call(.OpenMxLoaded, PACKAGE="OpenMx")
  mxSetDefaultOptions()
	if (.Platform$GUI!="Rgui") {
	  .Call(.enableMxLog)
	}
}

.onAttach <- function(libname, pkgname) {
	if (.Platform$GUI=="Rgui") {
		packageStartupMessage(paste("Notice: R GUI cannot display verbose output from the OpenMx backend.",
					    "If you need detail diagnostics then R CMD BATCH is one option."))
	}
	 if (.Platform$OS.type != "windows") {
		 if (!imxHasOpenMP()) {
			 packageStartupMessage("OpenMx may run faster if it is compiled to take advantage of multiple cores.")
		 } else if (Sys.getenv("OMP_NUM_THREADS") == "") {
			 packageStartupMessage(paste0("To take full advantage of multiple cores, use:\n",
				 "  mxOption(key='Number of Threads', value=parallel::detectCores()) #now\n",
				 "  Sys.setenv(OMP_NUM_THREADS=parallel::detectCores()) #before library(OpenMx)"))
		 }
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
#' The core function is \code{\link{mxModel}}, which makes a model. Models are containers for \code{\link{mxData}}, \code{\link[=mxMatrix]{matrices}}, \code{\link{mxPath}}s
#' \code{\link[=mxAlgebra]{algebras}}, \link{mxBounds}, \code{\link[=mxCI]{confidence intervals}}, and \code{\link[=mxConstraint]{mxConstraints}}.
#' Most models require an expectation (see the list below) to calculate the expectations for the model.
#' Models also need a fit function, several of which are built-in (see below).
#' OpenMx also allows user-defined fit functions for purposes not covered by the built-in functions. (e.g., \code{\link{mxFitFunctionR}} or \code{\link{mxFitFunctionAlgebra}}).
#' 
#' \emph{Note}, for mxModels of \code{type="RAM"}, the expectation and fit-function are set for you automatically.
#' 
#' \strong{Running and summarizing a model}
#' 
#' Once built, the resulting mxModel can be run (i.e., optimized) using  \code{\link{mxRun}}. This returns the fitted model.
#' 
#' You can summarize the results of the model using \code{\link[=summary.MxModel]{summary}}(yourModel)
#' 
#' \strong{Additional overview of model making and getting started}
#' 
#' The OpenMx manual is online (see references below). However, \code{\link{mxRun}}, \code{\link{mxModel}}, \code{\link{mxMatrix}}
#' all have working examples that will help get you started as well.
#'
#' The main OpenMx functions are: \code{\link{mxAlgebra}}, \code{\link{mxBounds}}, \code{\link{mxCI}}, \code{\link{mxConstraint}}, \code{\link{mxData}}, 
#' \code{\link{mxMatrix}}, \code{\link{mxModel}}, and \code{\link{mxPath}}.
#' 
#' Expectation functions include \code{\link{mxExpectationNormal}}, \code{\link{mxExpectationRAM}}, \code{\link{mxExpectationLISREL}}, and \code{\link{mxExpectationStateSpace}};
#' 
#' Fit functions include \link{mxFitFunctionML}, \link{mxFitFunctionAlgebra}, \link{mxFitFunctionRow} and \link{mxFitFunctionR}.
#' 
#' \strong{Datasets built into OpenMx}
#' 
#' \code{OpenMx} comes with over a dozen useful datasets built-in. Discover them using \code{data(package="OpenMx")}, and open them with, 
#' for example, \code{data("jointdata", package ="OpenMx", verbose= TRUE)}
#' 
#' Please cite the 'OpenMx' package in any publications that make use of it:
#' 
#' Michael C. Neale, Michael D. Hunter, Joshua N. Pritikin, Mahsa Zahery, Timothy R. Brick Robert M.
#' Kirkpatrick, Ryne Estabrook, Timothy C. Bates, Hermine H. Maes, Steven M. Boker. (2016).
#' OpenMx 2.0: Extended structural equation and statistical modeling. \emph{Psychometrika}, \strong{81}, 535–549. 
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
#' http://openmx.ssri.psu.edu/docs/OpenMx/latest/OpenMxUserGuide.pdf
#'
#' @references The OpenMx User's guide can be found at \url{https://openmx.ssri.psu.edu/documentation}
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
#' summary(m1)
#'
#' @docType package
#' @name OpenMx
NULL
