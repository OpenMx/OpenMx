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

removeTrailingSeparator <- function(x) {
	return(sub('/$', '', x))
}

#' Save model state to a checkpoint file
#'
#' @template args-model
#' @template args-chkpt
#'
#' @description
#' The function saves the last state of a model to a checkpoint file.
#' 
#' @details
#' In general, the arguments \sQuote{chkpt.directory} and \sQuote{chkpt.prefix} should be identical to the \code{\link{mxOption}}: \sQuote{Checkpoint Directory} and \sQuote{Checkpoint Prefix} that were specified on the model before execution.
#' 
#' Alternatively, the checkpoint file can be manually loaded as a data.frame in R.  Use \code{\link{read.table}} with the options \code{header=TRUE, sep="\t", stringsAsFactors=FALSE, check.names=FALSE}.
#'
#' @return
#' Returns a logical indicating the success of writing the checkpoint file to the checkpoint directory.
#' @template ref-manual
#' @family model state
#' @examples
#' library(OpenMx)
#' 
#' # Simulate some data
#' 
#' x=rnorm(1000, mean=0, sd=1)
#' y= 0.5*x + rnorm(1000, mean=0, sd=1)
#' tmpFrame <- data.frame(x, y)
#' tmpNames <- names(tmpFrame)
#' 
#' dir <- tempdir()  # safe place to create files
#' mxOption(key="Checkpoint Directory", value=dir)
#'
#' # Create a model that includes an expected covariance matrix,
#' # an expectation function, a fit function, and an observed covariance matrix
#' 
#' data <- mxData(cov(tmpFrame), type="cov", numObs = 1000)
#' expCov <- mxMatrix(type="Symm", nrow=2, ncol=2, values=c(.2,.1,.2), free=TRUE, name="expCov")
#' expFunction <- mxExpectationNormal(covariance="expCov", dimnames=tmpNames)
#' fitFunction <- mxFitFunctionML()
#' testModel <- mxModel(model="testModel", expCov, data, expFunction, fitFunction)
#' 
#' #Use mxRun to optimize the free parameters in the expected covariance matrix
#' modelOut <- mxRun(testModel)
#' modelOut$expCov
#' 
#' # Save the ending state of modelOut in a checkpoint file
#' mxSave(modelOut)
#' 
#' # Restore the saved model from the checkpoint file
#' modelSaved <- mxRestore(testModel)
#' modelSaved$expCov
mxSave <- function(model, chkpt.directory = ".", chkpt.prefix = "") {
	if (!is(model, "MxModel")) {
		stop("'model' argument must be a MxModel object")
	}
	if (!missing(chkpt.directory)) model <- mxOption(model,"Checkpoint Directory", chkpt.directory)
	if (!missing(chkpt.prefix))    model <- mxOption(model,"Checkpoint Prefix", chkpt.prefix)
	model <- mxOption(model,"Checkpoint Units",'evaluations')
	model <- mxOption(model,"Checkpoint Count",1)
	model <- mxModel(model, mxComputeOnce('fitfunction', 'fit'))
	mxRun(model, checkpoint=TRUE, silent=TRUE)
	invisible(TRUE)
}

#' Restore model state from a checkpoint file
#'
#' @template args-model
#' @template args-chkpt
#' @param line integer. Which line from the checkpoint file to restore (defaults to the last line)
#' @param strict logical. Require that the checkpoint name and model name match
#'
#' @details
#' In general, the arguments \sQuote{chkpt.directory} and \sQuote{chkpt.prefix} should be identical to the \code{\link{mxOption}}: \sQuote{Checkpoint Directory} and \sQuote{Checkpoint Prefix} that were specified on the model before execution.
#'
#' Alternatively, the checkpoint file can be manually loaded as a data.frame in R and passed to \code{\link{mxRestoreFromDataFrame}}.
#' Use \code{\link{read.table}} with the options \code{header=TRUE, sep="\t", stringsAsFactors=FALSE, check.names=FALSE}.
#'
#' @return
#' Returns an MxModel object with free parameters updated to the last
#' saved values. When \sQuote{line} is provided, the MxModel is updated
#' to the values on that line within the checkpoint file.
#'
#' @template ref-manual
#' @family model state
#' @examples
#' library(OpenMx)
#' 
#' # Simulate some data
#' 
#' x=rnorm(1000, mean=0, sd=1)
#' y= 0.5*x + rnorm(1000, mean=0, sd=1)
#' tmpFrame <- data.frame(x, y)
#' tmpNames <- names(tmpFrame)
#' 
#' dir <- tempdir()  # safe place to create files
#' mxOption(key="Checkpoint Directory", value=dir)
#' 
#' # Create a model that includes an expected covariance matrix,
#' # an expectation function, a fit function, and an observed covariance matrix
#' 
#' data <- mxData(cov(tmpFrame), type="cov", numObs = 1000)
#' expCov <- mxMatrix(type="Symm", nrow=2, ncol=2, values=c(.2,.1,.2), free=TRUE, name="expCov")
#' expFunction <- mxExpectationNormal(covariance="expCov", dimnames=tmpNames)
#' fitFunction <- mxFitFunctionML()
#' testModel <- mxModel(model="testModel", expCov, data, expFunction, fitFunction)
#' 
#' #Use mxRun to optimize the free parameters in the expected covariance matrix
#' modelOut <- mxRun(testModel, checkpoint = TRUE)
#' modelOut$expCov
#' 
#' #Use mxRestore to load the last checkpoint saved state of the model
#' modelRestore <- mxRestore(testModel)
#' modelRestore$expCov

mxRestore <- function(model, chkpt.directory = mxOption(model, "Checkpoint directory"),
		      chkpt.prefix = mxOption(model,"Checkpoint Prefix"), line=NULL, strict=FALSE)
{
  warnModelCreatedByOldVersion(model)
	chkpt.directory <- removeTrailingSeparator(chkpt.directory)
	if (strict && chkpt.prefix == "") {
		chkpt.prefix <- model$name
	}
	pattern <- paste("^\\Q", chkpt.prefix, "\\E.*(\\.omx)$", sep = '')
	chkpt.files <- list.files(chkpt.directory, full.names = FALSE)
	chkpt.files <- grep(pattern, chkpt.files, perl=TRUE, value=TRUE)
	if (strict) {
		if (length(chkpt.files) > 1) {
			stop(paste("chkpt.prefix", omxQuotes(chkpt.prefix),
				   "matched more than one file"))
		} else if (length(chkpt.files) == 0) {
			stop(paste("Cannot find", omxQuotes(paste(model$name, 'omx', sep=".")),
				   "in", chkpt.directory))
		}
	} else {
		if(length(chkpt.files) == 0) {
			return(model)
		}
		# Move the most likely match to the end so those estimates take precedence.
		matchIndex <- match(paste(model$name, 'omx', sep="."), chkpt.files)
		if (!is.na(matchIndex)) {
			chkpt.files <- c(chkpt.files[-matchIndex], paste(model$name, 'omx', sep="."))
		}
	}
	if (length(chkpt.files) > 1 && !is.null(line)) {
		stop(paste("Ambiguous: cannot specify line =", line,
			   "with more than one checkpoint found:",
			   omxQuotes(chkpt.files)))
	}
	if (length(chkpt.files) > 1) {
		message(paste("Loading estimates from more than one checkpoint:",
			      omxQuotes(chkpt.files)))
	}
	for(i in 1:length(chkpt.files)) {
		filename <- chkpt.files[[i]]
		filepath <- paste(chkpt.directory, filename, sep = '/')
		checkpoint <- read.table(filepath, header=TRUE, stringsAsFactors=FALSE, check.names=FALSE, sep="\t")
		model <- mxRestoreFromDataFrame(model, checkpoint, line)
	}
	return(model)
}

#' @rdname mxRestore
#' @param checkpoint a data.frame containing the model state
mxRestoreFromDataFrame <- function(model, checkpoint, line=NULL)
{
  warnModelCreatedByOldVersion(model)
  allPar <- names(omxGetParameters(model, indep=TRUE, free=NA))
  badPar <- intersect(c("OpenMxContext","OpenMxNumFree","OpenMxEvals",
			"iterations","timestamp","objective"), allPar)
  if (any(badPar)) {
    warning(paste("You model", omxQuotes(model$name), "has illegal parameters",
		  "names:", omxQuotes(badPar)))
  }
  pick <- match(allPar, colnames(checkpoint))
  if (all(is.na(pick))) return(model)
  allVars <- !any(is.na(pick))
  if (is.null(line)) {
    row <- nrow(checkpoint)
  } else {
    row <- line
    if (row < 2 || row > nrow(checkpoint)) {
      warning(paste("Requested line", line,
		    "but checkpoint contains lines 2 to",
		    nrow(checkpoint), "; using the last line"))
      row <- nrow(checkpoint)
    }
  }
  pick <- pick[!is.na(pick)]
  varNames <- colnames(checkpoint)[pick]
  model <- omxSetParameters(model, labels=varNames, values=as.numeric(checkpoint[row, varNames]))
  if (allVars) {
    model@output$fit <- checkpoint[row, 'objective']
    model@output$fitUnits <- checkpoint[row, 'fitUnits']
    if (!is.null(model@output$fitUnits) && model@output$fitUnits == '-2lnL') {
      model@output$Minus2LogLikelihood <- model@output$fit
    }
    model@output$estimates <- coef(model)
    se <- paste0(varNames, "SE")
    if (all(!is.na(match(se, colnames(checkpoint))))) {
      model@output$standardErrors <- matrix(as.numeric(checkpoint[row, se]), ncol=1,
					    dimnames=list(varNames,NULL))
    }
    grad <- paste0(varNames, "Grad")
    if (all(!is.na(match(grad, colnames(checkpoint))))) {
      vec <- as.numeric(checkpoint[row, grad])
      names(vec) <- varNames
      model@output$gradient <- vec
    }
    vc <- rep('', length(varNames) * (length(varNames)+1) / 2)
    vx <- 1
    for (cx in 1:length(varNames)) {
      for (rx in cx:length(varNames)) {
	vc[vx] <- paste0('V', varNames[rx], ':', varNames[cx])
	vx <- vx + 1
      }
    }
    if (all(!is.na(match(vc, colnames(checkpoint))))) {
      vcov <- matrix(0, length(varNames), length(varNames),
		     dimnames=list(varNames, varNames))
      vcov[lower.tri(vcov, TRUE)] <- as.numeric(checkpoint[row, vc])
      vcov[upper.tri(vcov)] <- t(vcov)[upper.tri(vcov)]
      model@output$vcov <- vcov
    }

    # refer to MxRun.R
    unsafe <- TRUE
    namespace <- imxGenerateNamespace(model)
    flatModel <- imxFlattenModel(model, namespace, unsafe)
    dependencies <- cycleDetection(flatModel)
    dependencies <- transitiveClosure(flatModel, dependencies)
    freeVarGroups <- buildFreeVarGroupList(flatModel)
    flatModel <- generateParameterList(flatModel, dependencies, freeVarGroups)
    matrices <- generateMatrixList(flatModel)
    parameters <- flatModel@parameters

    runstate <- model@runstate
    runstate$parameters <- parameters
    runstate$matrices <- matrices
    model@runstate <- runstate
    model@.wasRun <- TRUE
    model@.modifiedSinceRun <- FALSE #?
  }
  model
}
